#include "gbCnf.h"
#include "KNHome.h"

void KNHome::getBondChange(gbCnf& cnfModifier) {

  vector<double> direction = vdparams["direction"]; // e.g. 0.5 0.5 0
  vector<int> plane = viparams["plane"];
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
                                      bondEnergy, \
                                      oldBond, \
                                      cfg, \
                                      rng, \
                                      gen, \
                                      elems));
  }

  // do some statistics here

  // end do some statistics here

}

double KNHome::getBondChangeSingle(gbCnf& cnfModifier, \
                      const vector<double>& direction, \
                      const vector<int>& plane,
                      const double& LC, \
                      unordered_map<string, double>& bondEnergy, \
                      unordered_map<string, int>& oldBond, \
                      Config& cfg, \
                      std::default_random_engine& rng, \
                      std::uniform_int_distribution<int>& gen, \
                      const vector<string>& elems) {

  // randomly choose an atom
  int atomID = gen(rng);
  // shift to center, wrap, find vector, shift agian, wrap
  // assign sign on plane and direction randomly, also shuffle direction
  Config& cfgNew = cfg;
  cnfModifier.shiftAtomToCentral(cfgNew, atomID);
  auto newPlane = plane;
  auto newDirection = direction;
  for (auto& number:newPlane) {
    number *= (gen(rng) % 2 == 0) ? 1 : -1;
  }
  for (auto& number:newDirection) {
    number *= (gen(rng) % 2 == 0) ? 1 : -1;
  }
  unordered_map<string, int> newBond = cnfModifier.bondCountAll(cfgNew);
  // 2.860954037
  // 2.335959189
  // calculate difference in bonds
  double res = 0.0;
  for (auto&& i : newBond) {
    const string key = i.first;
    res += static_cast<double>(i.second - oldBond[key]) * bondEnergy[key];
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