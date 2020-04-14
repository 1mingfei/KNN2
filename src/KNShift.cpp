#include "gbCnf.h"
#include "KNHome.h"
#include "gbUtl.h"
void KNHome::getBondChange(gbCnf& cnfModifier) {

  vector<double> direction = vdparams["direction"]; // e.g. 0.5 0.5 0
  vector<double> plane = vdparams["plane"];
  double LC = dparams["LC"];
  vector<string> bondName = vsparams["bondName"];
  vector<double> bondEnergyInit = vdparams["bondEnergy"];
  unordered_map<string, double> bondEnergy;
  // bond Energy
  for (int i = 0; i < bondName.size(); ++i) {
    bondEnergy[bondName[i]] = bondEnergyInit[i];
  }

  string fname = sparams["initconfig"];
  Config cfg = std::move(cnfModifier.readCfg(fname));
  int N = iparams["NBarriers"];
  vector<string> elems = vsparams["elems"];
  std::default_random_engine rng(std::random_device{}() + me);
  std::uniform_int_distribution<int> gen(0, cfg.natoms - 1);
  vector<double> res;
  vector<int> factors = viparams["factors"];

  MPI_Barrier(MPI_COMM_WORLD);
  if (nProcs == 1)
    cnfModifier.getNBL_serial(cfg, 3.5);
  else
    cnfModifier.getNBL(cfg, 3.5);

  unordered_map<string, int> oldBond = cnfModifier.bondCountAll(cfg);

  for (int i = 0; i < N; ++i) {
    res.push_back(getBondChangeSingle(cnfModifier, \
                                      direction, \
                                      plane, \
                                      LC, \
                                      factors,\
                                      bondEnergy, \
                                      cfg, \
                                      rng, \
                                      gen));
  }

  // do some statistics here

  // end do some statistics here

}

double KNHome::getBondChangeSingle(gbCnf& cnfModifier, \
                      const vector<double>& direction, \
                      const vector<double>& plane,
                      const double& LC, \
                      const vector<int>& factors, \
                      unordered_map<string, double>& bondEnergy, \
                      Config& cfg, \
                      std::default_random_engine& rng, \
                      std::uniform_int_distribution<int>& gen) {

  // randomly choose an atom
  int atomID = gen(rng);
  // shift to center, wrap, find vector, shift agian, wrap
  // assign sign on plane and direction randomly, also shuffle direction
  Config& cfgNew = cfg;
  cnfModifier.shiftAtomToCentral(cfgNew, atomID);
  auto newPlane = plane;
  for (auto& number:newPlane) {
    number *= (gen(rng) % 2 == 0) ? 1 : -1;
  }
  vector<double> planeMoveDistance;
  planeMoveDistance.reserve(3);
  double newPlaneInner =  vecInnProd33(newPlane, newPlane);
  for (int i = 0; i < 3; ++i) {
    // 2 here means only move a half distance to make sure atoms not aross
    // the center
    planeMoveDistance.push_back(
        newPlane[i] * LC / newPlaneInner / 2);
  }
  cnfModifier.shiftPst(cfgNew, planeMoveDistance);

  // Calculate the planes, there should be three planes, move half of the atoms
  // according to the first plane, calculate bonds change between second and
  // third plane.
  // First should be Ax + By + Cz = D = Lx/2 + Ly/2 + Lx/2. Second and third
  // should be Ax + By + Cz = D = Lx/2 + Ly/2 + Lx/2 ± LC/sqrt(h^2+k^2+l^2)

  double D1 = 0.5*(cfgNew.bvx[0] + cfgNew.bvy[1] + cfgNew.bvy[2]);
  double D2 = D1 + LC/sqrt(newPlaneInner);
  double D3 = D1 - LC/sqrt(newPlaneInner);

  double checkResult;
  // Store atoms' index between between second plane and third plane.
  vector<int> nearAtomList;
  for (const auto& atom:cfgNew.atoms) {
    checkResult = atom.pst[0] * newPlane[0] + atom.pst[1] * newPlane[1]
        + atom.pst[2] * newPlane[2];
    if (checkResult > D2 || checkResult < D3)
      continue;
    nearAtomList.push_back(atom.id);
  }
  // Calculate original bonds number
  unordered_map<string, int> bondsCountBefore;
  string tp1, tp2, bond;
  for (auto it1 = nearAtomList.begin(); it1 < nearAtomList.end(); ++it1) {
    tp1 = cfgNew.atoms[*it1].tp;
    for (auto it2 = nearAtomList.begin(); it2 < it1; ++it2) {
      tp2 = cfgNew.atoms[*it2].tp;

      if (tp1.compare(tp2) < 0) {
        bond = tp1;
        bond += "-";
        bond += tp2;
      } else {
        bond = tp2;
        bond += "-";
        bond += tp1;
      }
      if (cnfModifier.calDist({cfgNew.bvx[0], cfgNew.bvy[1], cfgNew.bvy[2]},
                              cfgNew.atoms[*it1],
                              cfgNew.atoms[*it2]) < 3.2)
        bondsCountBefore[bond]++;
    }
  }
  // move atoms
  auto newDirection = direction;
  for (auto& number:newDirection) {
    number *= (gen(rng) % 2 == 0) ? 1 : -1;
  }
  vector<double> atomMoveDistance;
  atomMoveDistance.reserve(3);
  for (int i = 0; i < 3; ++i) {
    atomMoveDistance.push_back(LC*newDirection[i]);
  }
  for(const auto &index:nearAtomList){
    checkResult = cfgNew.atoms[index].pst[0] * newPlane[0]
        + cfgNew.atoms[index].pst[1] * newPlane[1]
        + cfgNew.atoms[index].pst[2] * newPlane[2];
    if (checkResult < D1)
      continue;
    for (int i = 0; i < 3; ++i) {
      cfgNew.atoms[index].pst[i] += atomMoveDistance[i];
    }
  }
  cnfModifier.wrapAtomPos(cfgNew);

  // Calculate bonds number after moving
  unordered_map<string, int> bondsCountAfter;
  for (auto it1 = nearAtomList.begin(); it1 < nearAtomList.end(); ++it1) {
    tp1 = cfgNew.atoms[*it1].tp;
    for (auto it2 = nearAtomList.begin(); it2 < it1; ++it2) {
      tp2 = cfgNew.atoms[*it2].tp;
      if (tp1.compare(tp2) < 0) {
        bond = tp1;
        bond += "-";
        bond += tp2;
      } else {
        bond = tp2;
        bond += "-";
        bond += tp1;
      }
      if (cnfModifier.calDist({cfgNew.bvx[0], cfgNew.bvy[1], cfgNew.bvy[2]},
                              cfgNew.atoms[*it1],
                              cfgNew.atoms[*it2]) < 3.2)
        bondsCountAfter[bond]++;
    }
  }

  // unordered_map<string, int> newBond = cnfModifier.bondCountAll(cfgNew);

  // calculate difference in bonds
  double res = 0.0;
  for (const auto& [key, count] : bondsCountAfter) {
    res += static_cast<double>(count - bondsCountBefore[key]) * bondEnergy[key];
  }
  return res;
}

unordered_map<string, int> gbCnf::bondCountAll(const Config& cfg) {
  unordered_map<string, int> bondsCount;
  string tp1, tp2, bond;
  for (const auto& atom:cfg.atoms) {
    tp1 = atom.tp;
    for (const auto& atom2ID:atom.FNNL) {
      tp2 = cfg.atoms[atom2ID].tp;
      if (tp1.compare(tp2) < 0) {
        bond = tp1;
        bond += "-";
        bond += tp2;
      } else {
        bond = tp2;
        bond += "-";
        bond += tp1;
      }
      bondsCount[bond]++;
    }
  }
  for (auto &[bond, count]:bondsCount) {
    count /= 2;
  }
  return bondsCount;
}

void gbCnf::shiftAtomToCentral(Config& cnf, const int& id) {
  array<double, 3> criterionDistance =
      {cnf.atoms[id].prl[0], cnf.atoms[id].prl[1], cnf.atoms[id].prl[2]};
  for (auto& atm:cnf.atoms) {
    for (int i = 0; i < 3; ++i) {
      atm.prl[i] -= criterionDistance[i];
    }
  }
  wrapAtomPrl(cnf);
  cnvprl2pst(cnf);
}
void gbCnf::shiftPst(Config& cnf, const vector<double>& distance){
  for (auto& atm:cnf.atoms) {
    for (int i = 0; i < 3; ++i) {
      atm.pst[i] -= distance[i];
    }
  }
  wrapAtomPos(cnf);
  cnvpst2prl(cnf);
}