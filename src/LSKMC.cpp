#include "gbCnf.h"
#include "KNHome.h"
#include "LSKMC.h"
#include "KNUtility.h"

#define KB 8.6173303e-5
#define KB_INV 11604.5221105
#define NEI_NUMBER 12
#define EPSILON 1e-15

namespace LS {

LSKMC::LSKMC(gbCnf& cnfModifierIn, \
             Config& c0In, \
             unordered_map<string, double>& embeddingIn, \
             vector<int>& vacListIn, \
             string& EDiffIn, \
             Model& k2pModelBIn, \
             Model& k2pModelDIn, \
             double& RCutIn, \
             double& RCut2In, \
             double& temperatureIn, \
             double& timeIn, \
             double& prefixIn, \
             double& ECutoffIn, \
             long long& stepIn, \
             ofstream& ofsIn, \
             bool& switchLSKMCIn, \
             bool& switchUnknownIn, \
             vector<string>& elemsIn, \
             vector<double>& elemsEffectOffsetIn)
  : cnfModifier(cnfModifierIn), \
    c0(c0In), \
    embedding(embeddingIn), \
    vacList(vacListIn), \
    EDiff(EDiffIn), \
    k2pModelB(k2pModelBIn), \
    k2pModelD(k2pModelDIn), \
    RCut(RCutIn), \
    RCut2(RCut2In), \
    temperature(temperatureIn), \
    time(timeIn), \
    prefix(prefixIn), \
    ECutoff(ECutoffIn), \
    step(stepIn), \
    ofs(ofsIn), \
    switchLSKMC(switchLSKMCIn), \
    switchUnknown(switchUnknownIn), \
    elems(elemsIn), \
    elemsEffectOffset(elemsEffectOffsetIn)
{

  eventMap.clear();
  // watch out
  // this function need to be updated if multiple vacacies is in use
  searchStatesDFS();

}

void LSKMC::testCnfModification() {

  cnfModifier.writeCfgData(c0, "testOut.1.cfg");
  c0.atoms[0].tp = "X";
  cnfModifier.writeCfgData(c0, "testOut.2.cfg");

}

void LSKMC::test_vvd2mat() {

  vvd vIn = { {0.1234, 0.234}, \
              {0.0223, 0.388}, \
              {234.0, 2134} };
  for (int i = 0; i < vIn.size(); ++i) {
    for (int j = 0; j < vIn[0].size(); ++j) {
      cout << vIn[i][j] << " ";
    }
    cout << endl;
  }
  cout << endl << endl;
  mat Arm = vvd2mat(vIn);
  Arm.print();

}

void LSKMC::helperDFS(const int& i, \
               const int& vac, \
               unordered_set<int>& visited) {

  if (visited.count(i))
    return;
  visited.insert(i);

  for (auto&& j : c0.atoms[i].FNNL) {
    if (visited.count(j))
      continue;

    getOrPutEvent(i, j);
    string tmpHash = to_string(i) + "_" + to_string(j);
    double currBarr = eventMap[tmpHash].getBarrier();
    barriers.push_back(currBarr);
    if (currBarr < ECutoff) {
      trapList[vac].insert(j);
      helperDFS(j, vac, visited);
    } else {
      absorbList[vac].insert(j);
    }
  }

  return;
}

void LSKMC::searchStatesDFS() {

  for (int i = 0; i < vacList.size(); ++i) {
    unordered_set<int> visited;
    helperDFS(vacList[i], vacList[i], visited);

#ifdef DEBUG_TRAP
    cout << "vac # " << vacList[i] << " trap size : " \
         << trapList[vacList[i]].size() << endl;
    cout << "vac # " << vacList[i] << " absortb size : " \
         << absorbList[vacList[i]].size() << endl;
#endif
  }
  barrierStats();
}

void LSKMC::outputTrapCfg(const int& vac, const string& fname) {

  Config cNew;
  cNew.length = c0.length;
  cNew.cell = c0.cell;
  // bvx, tvx, bvy, tvy, bvz, tvz;
  cNew.bvx = c0.bvx;
  cNew.tvx = c0.tvx;
  cNew.bvy = c0.bvy;
  cNew.tvy = c0.tvy;
  cNew.bvz = c0.bvz;
  cNew.tvz = c0.tvz;
  cNew.natoms = trapList[vac].size();
  for (const auto& i : trapList[vac]) {
    cNew.atoms.push_back(c0.atoms[i]);
  }
  cnfModifier.writeCfgData(cNew, fname);

}

void LSKMC::outputAbsorbCfg(const int& vac, const string& fname) {

  Config cNew;
  cNew.length = c0.length;
  cNew.cell = c0.cell;
  // bvx, tvx, bvy, tvy, bvz, tvz;
  cNew.bvx = c0.bvx;
  cNew.tvx = c0.tvx;
  cNew.bvy = c0.bvy;
  cNew.tvy = c0.tvy;
  cNew.bvz = c0.bvz;
  cNew.tvz = c0.tvz;
  cNew.natoms = absorbList[vac].size();
  for (const auto& i : absorbList[vac]) {
    cNew.atoms.push_back(c0.atoms[i]);
  }
  cnfModifier.writeCfgData(cNew, fname);

}

void LSKMC::barrierStats() {
  double sum = std::accumulate(barriers.begin(), barriers.end(), 0.0);
  std::sort(barriers.begin(), barriers.end());
  ofstream ofs("barrier_stats.txt", std::ofstream::out | std::ofstream::app);
  ofs << "step " << step \
      << " mean: " << (sum / static_cast<double>(barriers.size())) \
      << " min: " << barriers[0] \
      << " 25%: " << barriers[static_cast<int>(barriers.size() / 4.0) - 1] \
      << " 50%: " << barriers[static_cast<int>(barriers.size() / 2.0) - 1] \
      << " 75%: " << barriers[static_cast<int>(barriers.size() * 0.75) - 1] \
      << " max: " << barriers[barriers.size() - 1] \
      << endl;
  ofs.close();
}

void LSKMC::calVVD_M(const int& vac) {

  int size = trapList[vac].size() + absorbList[vac].size();
  VVD_M.resize(size);
  for (int i = 0; i < size; ++i)
    VVD_M[i].resize(size, 0.0);

  // build maps forward and backward
  mapAtomID2MatID.clear();
  mapMatID2AtomID.clear();
  int id = 0;

#ifdef DEBUG_TRAP
  cout << "matrix size : " << size << endl;
#endif

  // absorption states MUST come ahead of transient states
  for (const auto& i : absorbList[vac]) {
    mapAtomID2MatID[i] = id;
    mapMatID2AtomID[id] = i;
    VVD_M[id][id] = 1.0;
    ++id;
  }

  for (const auto& i : trapList[vac]) {
    mapAtomID2MatID[i] = id;
    mapMatID2AtomID[id] = i;
    VVD_M[id][id] = 0.0;
    ++id;
  }

  // calculate Tau vector
  getVD_Tau(vac, size);

  // take care of i != j
  // i, j: matrix row and col
  // id_i, id_j: atom id
  for (int i = absorbList[vac].size(); i < size; ++i) {
    int id_i = mapMatID2AtomID[i];
    for (int j = 0; j < size; ++j) {
      if (i == j)
        continue;

      int id_j = mapMatID2AtomID[j];
      if (find(c0.atoms[id_i].FNNL.begin(), c0.atoms[id_i].FNNL.end(), id_j) \
          != c0.atoms[id_i].FNNL.end()) {
        getOrPutEvent(id_i, id_j);
        string tmpHash = to_string(id_i) + "_" + to_string(id_j);
        VVD_M[i][j] = eventMap[tmpHash].getRate() \
                      * VD_Tau[i - absorbList[vac].size()]; // correct
      }
    }
  }

  Arm_M = vvd2mat(VVD_M);

#ifdef DEBUG_TRAP
  cout << "M : " << VVD_M.size() << "x" << VVD_M[0].size() << endl;
  // Arm_M.print();
  cout.width(4);
  cout << Arm_M << endl;
#endif

}

void LSKMC::getVD_Tau(const int& vac, const int& size) {

  // i, j: matrix row and col
  // id_i, id_j: atom id

  VD_Tau.clear();
  VD_Tau.resize(trapList[vac].size(), 0.0);
  for (int i = absorbList[vac].size(); i < size; ++i) {
    double res = 0.0;
    int id_i = mapMatID2AtomID[i];
    for (const int& id_j : c0.atoms[id_i].FNNL) {
      getOrPutEvent(id_i, id_j);
      string tmpHash = to_string(id_i) + "_" + to_string(id_j);
      res += eventMap[tmpHash].getRate();
    }

    if (res < EPSILON)
      res += EPSILON;
    VD_Tau[i - absorbList[vac].size()] = 1.0 / res;
  }
  Arm_Tau = vd2vec(VD_Tau);

#ifdef DEBUG_TRAP
  cout << "Tau : " << VD_Tau.size() << endl;
  cout.width(4);
  cout << Arm_Tau << endl;
#endif

}

void LSKMC::calVVD_R(const int& vac) {
  VVD_R.resize(trapList[vac].size());
  for (int i = absorbList[vac].size(); i < VVD_M.size(); ++i) {
    vector<double> tmp(absorbList[vac].size(), 0.0);
    for (int j = 0; j < absorbList[vac].size(); ++j) {
      tmp[j] = VVD_M[i][j];
    }
    VVD_R[i - absorbList[vac].size()] = tmp;
  }
  Arm_R = vvd2mat(VVD_R);

#ifdef DEBUG_TRAP
  cout << "R : " << VVD_R.size() << "x" << VVD_R[0].size() << endl;
  cout.width(4);
  cout << Arm_R << endl;
#endif

}

void LSKMC::calVVD_T(const int& vac) {
  VVD_T.resize(trapList[vac].size());
  for (int i = absorbList[vac].size(); i < VVD_M.size(); ++i) {
    vector<double> tmp(trapList[vac].size(), 0.0);
    for (int j = absorbList[vac].size(); j < VVD_M.size(); ++j) {
      tmp[j - absorbList[vac].size()] = VVD_M[i][j];
    }
    VVD_T[i - absorbList[vac].size()] = tmp;
  }
  Arm_T = vvd2mat(VVD_T);

#ifdef DEBUG_TRAP
  cout << "T : " << VVD_T.size() << "x" << VVD_T[0].size() << endl;
  cout.width(4);
  cout << Arm_T << endl;
#endif

}

void LSKMC::getOrPutEvent(const int& i, const int& j) {
  string tmpHash = to_string(i) + "_" + to_string(j);
  if (eventMap.find(tmpHash) == eventMap.end()) {
    vector<double> currBarr = cnfModifier.calBarrierAndEdiff(c0, \
                                          temperature, \
                                          RCut2, \
                                          EDiff, \
                                          embedding, \
                                          k2pModelB, \
                                          k2pModelD, \
                                          make_pair(i, j), \
                                          switchUnknown, \
                                          elems, \
                                          elemsEffectOffset);
    LSEvent event(make_pair(i, j));
    event.setBarrier(currBarr[0]);
    event.setRate(prefix * exp(-currBarr[0] * KB_INV / temperature));
    event.setEnergyChange(currBarr[1]);
    eventMap[tmpHash] = event;
  }
}

void LSKMC::calExitTimePi(const int& vac) {

  calVVD_M(vac);
  calVVD_R(vac);
  calVVD_T(vac);

  vd VD_P0(VVD_T.size(), 0.0);
  VD_P0[mapAtomID2MatID[vac]] = 1.0;
  vec Arma_P0 = vd2vec(VD_P0);

#ifdef DEBUG_TRAP
  cout << "Arma_P0 : " << endl << Arma_P0 << endl;
#endif

  mat sharedMatrix = Arma_P0.t() \
                     * arma::pinv(arma::eye(arma::size(Arm_T)) - Arm_T);
  exitTime = arma::as_scalar(sharedMatrix * Arm_Tau);

#ifdef DEBUG_TRAP
  cout << "exit time : " << exitTime << endl;
#endif

  Arm_Pi = sharedMatrix * Arm_R;
  double sum = accu(Arm_Pi);
  Arm_Pi /= sum;

#ifdef DEBUG_TRAP
  cout << "Pi : " << endl;
  cout << Arm_Pi << endl;
  cout << Arm_Pi.n_rows << "x" << Arm_Pi.n_cols << endl;
  sum = accu(Arm_Pi);
  cout << "total (should be 1.0) :" << sum << endl;
#endif
}

void LSKMC::updateTime() {

  if (exitTime > 0.0) {
    time += exitTime;
  }

}

bool LSKMC::validTrap(const int& vac) {
  if (trapList[vac].size() == 0)
    return false;
  if (absorbList[vac].size() == 0)
    return false;
  calExitTimePi(vac);
  if (exitTime > 1.0e5)
    return false;
  return true;
}

void LSKMC::selectAndExecute(const int& vac) {

  if (!validTrap(vac))
    return;

  double randVal = (double) rand() / (RAND_MAX);
  vd prob = mat2vd(Arm_Pi);
  vd probAccu(prob.size(), 0.0);
  for (int i = 0; i < prob.size(); ++i) {
    if (i == 0)
      probAccu[i] = prob[i];
    else {
      probAccu[i] = probAccu[i - 1] + prob[i];
    }
  }

  auto it = lower_bound(probAccu.begin(), probAccu.end(), randVal);

  int dist = 0;
  if (it == probAccu.cend()) {
    dist = probAccu.size() - 1;
  }

  dist = distance(probAccu.begin(), it);
  int iFirst = vac;
  // starting from absorbing state, no offset needed
  int iSecond = mapMatID2AtomID[dist];
  LSEvent lsevent(make_pair(iFirst, iSecond));

  if (exitTime > 4.0)
    cnfModifier.writeCfgData(c0, "long_lskmc_iter_" \
                                 + to_string(step) + "_0.cfg");

  lsevent.exeEvent(c0, RCut);

  if (exitTime > 4.0)
    cnfModifier.writeCfgData(c0, "long_lskmc_iter_" \
                                 + to_string(step) + "_1.cfg");

  updateTime();
  ofs << "# LSKMC " << step << " " << time << " ave exit time : " \
       << exitTime << endl;

  if (switchLSKMC) {
    cnfModifier.writeCfgData(c0, "lskmc_iter_" + to_string(step) + ".cfg");
  }

#ifdef DEBUG_SELECT_TRAP
  for (int i = 0; i < probAccu.size(); ++i)
    cout << i << " " << probAccu[i] << endl;
  cout << "random number : " << randVal << endl;
  cout << dist << endl;
#endif


}

} // end namespace LS


void KNHome::LSKMCOneRun(gbCnf& cnfModifier) {

  KMCInit(cnfModifier);
  buildEventList(cnfModifier);
  LS::LSKMC lskmc(cnfModifier, \
                  c0, \
                  embedding, \
                  vacList, \
                  EDiff, \
                  k2pModelB, \
                  k2pModelD, \
                  RCut, \
                  RCut2, \
                  temperature, \
                  time, \
                  prefix, \
                  ECutoff, \
                  step, \
                  ofs, \
                  switchLSKMC, \
                  switchUnknown, \
                  elems, \
                  elemsEffectOffset);

#ifdef DEBUG_TRAP
  for (const auto& i : vacList) {
    lskmc.outputAbsorbCfg(i, "debug_absorb_" + to_string(i) + ".cfg");
    lskmc.outputTrapCfg(i, "debug_trap_" + to_string(i) + ".cfg");
    lskmc.barrierStats();
    if (validTrap()) {
      lskmc.selectAndExecute(i);
    }
  }
#endif

}

bool KNHome::isTrapped(const double& oneStepTime) {
  long long trapStepCriteria = iparams["trapStepCriteria"];
  double oneTrapTimeCriteria = dparams["oneTrapTimeCriteria"];
  if (oneStepTime < oneTrapTimeCriteria)
    ++trapStep;
  else
    trapStep = 0;
  if (trapStep < trapStepCriteria)
    return false;
  else {
    trapStep = 0;
    return true;
  }
}

void KNHome::LSKMCSimulation(gbCnf& cnfModifier) {

  KMCInit(cnfModifier);

  double oneStepTime = 0.0;

  while (iter < maxIter) {

    if (nProcs == 1)
      buildEventList_serial(cnfModifier);
    else
      buildEventList(cnfModifier);

    int eventID = 0;
    if (me == 0) {
      auto&& event = selectEvent(eventID);
    }

    MPI_Bcast(&eventID, 1, MPI_INT, 0, MPI_COMM_WORLD);

    eventList[eventID].exeEvent(c0, RCut); // event updated

    if (me == 0) {
      double oneStepTime = updateTime();
      updateEnergy(eventID);
    }

    ++step;
    ++iter;

    if (me == 0) {
      if (isTrapped(oneStepTime)) {
        LS::LSKMC lskmc(cnfModifier, \
                        c0, \
                        embedding, \
                        vacList, \
                        EDiff, \
                        k2pModelB, \
                        k2pModelD, \
                        RCut, \
                        RCut2, \
                        temperature, \
                        time, \
                        prefix, \
                        ECutoff, \
                        step, \
                        ofs,\
                        switchLSKMC, \
                        switchUnknown, \
                        elems, \
                        elemsEffectOffset);
        for (const auto& i : vacList) {
          lskmc.selectAndExecute(i);
        }
      }

      if (step % nTallyOutput == 0)
      ofs << std::setprecision(7) << step << " " << time << " " \
          << E_tot << " " << endl;

      if (step % nTallyConf == 0)
        cnfModifier.writeCfgData(c0, to_string(step) + ".cfg");
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

}