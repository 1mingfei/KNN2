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
             Model& k2pModelBIn, \
             Model& k2pModelDIn, \
             double& RCutIn, \
             double& RCut2In, \
             double& temperatureIn, \
             double& timeIn, \
             double& prefixIn, \
             double& E_totIn, \
             double& ECutoffIn, \
             long long& maxIterIn, \
             long long& iterIn, \
             long long& stepIn, \
             int& nTallyConfIn, \
             int& nTallyOutput)
  : cnfModifier(cnfModifierIn), \
    c0(c0In), \
    embedding(embeddingIn), \
    vacList(vacListIn), \
    k2pModelB(k2pModelBIn), \
    k2pModelD(k2pModelDIn), \
    RCut(RCutIn), \
    RCut2(RCut2In), \
    temperature(temperatureIn), \
    time(timeIn), \
    prefix(prefixIn), \
    E_tot(E_totIn), \
    ECutoff(ECutoffIn), \
    maxIter(maxIterIn), \
    iter(iterIn), \
    step(stepIn), \
    nTallyConf(nTallyConfIn), \
    nTallyOutput(nTallyOutput)
{
  // initialize eventMap
  eventMap.clear();
}

void LSKMC::testCnfModification() {
  cnfModifier.writeCfgData(c0, "testOut.1.cfg");
  c0.atoms[0].tp = "X";
  cnfModifier.writeCfgData(c0, "testOut.2.cfg");
  return;
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

    // vector<double> currBarr = cnfModifier.calBarrierAndEdiff(c0, \
    //                                       temperature, \
    //                                       RCut2, \
    //                                       false, \
    //                                       embedding, \
    //                                       k2pModelB, \
    //                                       k2pModelD, \
    //                                       make_pair(i, j));

    // string tmpHash = to_string(i) + "_" + to_string(j);
    // KMCEvent event(make_pair(i, j));
    // event.setBarrier(currBarr[0]);
    // event.setRate(exp(-currBarr[0] * KB_INV / temperature));
    // event.setEnergyChange(currBarr[1]);
    // eventMap[tmpHash] = event;
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
    cout << "vac # " << vacList[i] << " : " << trapList[vacList[i]].size() \
         << endl;
    cout << "vac # " << vacList[i] << " : " << absorbList[vacList[i]].size() \
         << endl;
  }
  return;
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
  return;
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
  return;
}

void LSKMC::barrierStats() {
  double sum = std::accumulate(barriers.begin(), barriers.end(), 0.0);
  std::sort(barriers.begin(), barriers.end());
  ofstream ofs("barrier_stats.txt", std::ofstream::out);
  ofs << "mean: " << (sum / static_cast<double>(barriers.size())) \
      << " min: " << barriers[0] \
      << " 25%: " << barriers[static_cast<int>(barriers.size() / 4.0) - 1] \
      << " 50%: " << barriers[static_cast<int>(barriers.size() / 2.0) - 1] \
      << " 75%: " << barriers[static_cast<int>(barriers.size() * 0.75) - 1] \
      << " max: " << barriers[barriers.size() - 1] \
      << endl;
}

void LSKMC::calVVD_M(const int& vac) {

  int size = trapList.size() + absorbList.size();
  VVD_M.resize(size);
  for (int i = 0; i < size; ++i)
    VVD_M[i].resize(size, 0.0);

  // calculate Tau vector
  getVD_Tau(size);

  // build maps forward and backward
  mapAtomID2MatID.clear();
  mapMatID2AtomID.clear();
  int id = 0;

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


  // take care of i != j
  // i, j: matrix row and col
  // id_i, id_j: atom id
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      if (i == j)
        continue;

      int id_i = mapMatID2AtomID[i], id_j = mapMatID2AtomID[j];
      if (find(c0.atoms[id_i].FNNL.begin(), c0.atoms[id_i].FNNL.end(), id_j) \
          != c0.atoms[id_i].FNNL.end()) {
        getOrPutEvent(id_i, id_j);
        string tmpHash = to_string(id_i) + "_" + to_string(id_j);
        VVD_M[id_i][id_j] = eventMap[tmpHash].getRate() * VD_Tau[i]; //
      }
    }
  }
}

void LSKMC::getVD_Tau(const int& size) {

  // i, j: matrix row and col
  // id_i, id_j: atom id

  VD_Tau.clear();
  for (int i = 0; i < size; ++i) {
    double res = 0.0;
    for (const int& j : c0.atoms[i].FNNL) {
      int id_i = mapMatID2AtomID[i], id_j = mapMatID2AtomID[j];
      getOrPutEvent(id_i, id_j);
      string tmpHash = to_string(id_i) + "_" + to_string(id_j);
      res += eventMap[tmpHash].getRate();
    }

    if (res < EPSILON)
      res += EPSILON;
    VD_Tau.push_back(1.0 / res);
  }
}

void LSKMC::getOrPutEvent(const int& i, const int& j) {
  string tmpHash = to_string(i) + "_" + to_string(j);
  if (eventMap.find(tmpHash) == eventMap.end()) {
    vector<double> currBarr = cnfModifier.calBarrierAndEdiff(c0, \
                                          temperature, \
                                          RCut2, \
                                          false, \
                                          embedding, \
                                          k2pModelB, \
                                          k2pModelD, \
                                          make_pair(i, j));
    KMCEvent event(make_pair(i, j));
    event.setBarrier(currBarr[0]);
    event.setRate(exp(-currBarr[0] * KB_INV / temperature));
    event.setEnergyChange(currBarr[1]);
    eventMap[tmpHash] = event;
  }
}

} // end namespace LS


void KNHome::LSKMCSimulation(gbCnf& cnfModifier) {
  KMCInit(cnfModifier);
  buildEventList(cnfModifier);
  LS::LSKMC lskmc(cnfModifier, \
                  c0, \
                  embedding, \
                  vacList, \
                  k2pModelB, \
                  k2pModelD, \
                  RCut, \
                  RCut2, \
                  temperature, \
                  time, \
                  prefix, \
                  E_tot, \
                  ECutoff, \
                  maxIter, \
                  iter, \
                  step, \
                  nTallyConf, \
                  nTallyOutput);

  lskmc.searchStatesDFS();

#ifdef DEBUG_TRAP
  for (const auto& i : vacList) {
    cout << "outputing : " << i << endl;
    lskmc.outputAbsorbCfg(i, "debug_absorb_" + to_string(i) + ".cfg");
    lskmc.outputTrapCfg(i, "debug_trap_" + to_string(i) + ".cfg");
    lskmc.barrierStats();
  }
#endif

}
