/*
 * Author: 1mingfei 
 * Date:   2019-10-23
 * Purpose: encoding based on Bond Count
 * self-explained
 */

#include "gbCnf.h"
#include "KNHome.h"
// #include "armadillo"
// using arma::mat;
// using arma::vec;



/* output hashmap to file */
template<class T, class Y>
inline void writeHashmap(const string& fname, const int i, const int j, \
                        const unordered_map<T, Y>& v, \
                        const bool isOutputKey = false) {
  ofstream ofs(fname, std::ofstream::app);
  ofs << "config " << i << " end " << j << " ";
  if (isOutputKey) {
    for (const auto& val : v) 
      ofs << std::setw(4) << val.first << " ";
  } else {
    for (const auto& val : v) 
      ofs << std::setw(4) << val.second << " ";
  }
  ofs << endl;
}

void KNHome::KNBondCount() {
  gbCnf cnfModifier(*this);
  vector<string> elems = vsparams["elems"];
  int NConfigs = iparams["NConfigs"];
  int NBars = iparams["NBarriers"];
  double RCut = dparams["RCut"];
  //string fpname = "pairs.txt";//
  string fpname = sparams["PairFile"];
  vector<vector<int>> pairs = readPairs(fpname);

  unordered_map<string, int> res;
  res["Al-Al"] = 0;
  res["Al-Mg"] = 0;
  res["Al-Zn"] = 0;
  res["Al-X"] = 0;
  res["Mg-Mg"] = 0;
  res["Mg-Zn"] = 0;
  res["Mg-X"] = 0;
  res["Zn-Zn"] = 0;  
  res["Zn-X"] = 0;
  res["X-X"] = 0;
  ofstream ofs("bondCount.txt", std::ofstream::app);
  for (const auto& val : res) 
    ofs << std::setw(4) << val.first << " ";
  ofs << "\n";

  MPI_Barrier(MPI_COMM_WORLD);
  if (me == 0)
    std::cout << "count bondings of the configurations\n";
  int quotient = pairs.size() / nProcs;
  int remainder = pairs.size() % nProcs;
  int nCycle = remainder ? (quotient + 1) : quotient;

  for (int j = 0; j < nCycle; ++j) {
    for (int i = (j * nProcs); i < ((j + 1) * nProcs); ++i) {
      if ((i % nProcs != me) || (i >= pairs.size())) continue;
      //string fname = "in.cfg";//
      string fname = "config" + to_string(pairs[i][0]) + "/s/start.cfg";

      Config cfg = cnfModifier.readCfg(fname);
      vector<int> pair = {pairs[i][2], pairs[i][3]};
      cnfModifier.getNBL(cfg, RCut);
      unordered_map<string, int> pairCount;
      pairCount = cnfModifier.countPairs(cfg, elems, pair);
      writeHashmap("bondCount.txt", pairs[i][0], pairs[i][1], pairCount, false);
    }
  }
}

unordered_map<string, int> KNHome::gbCnf::countPairs(Config& cnf, \
                                                 const vector<string>& elems, \
                                                 const vector<int>& pair) {
  unordered_map<string, int> res;
  res["Al-Al"] = 0;
  res["Al-Mg"] = 0;
  res["Al-Zn"] = 0;
  res["Al-X"] = 0;
  res["Mg-Mg"] = 0;
  res["Mg-Zn"] = 0;
  res["Mg-X"] = 0;
  res["Zn-Zn"] = 0;  
  res["Zn-X"] = 0;
  res["X-X"] = 0;
  for (int i = 0; i < 2; ++i) {
    KNAtom atm = cnf.atoms[pair[i]];
    for (const int ii : atm.NBL) {
      KNAtom nbAtm = cnf.atoms[ii];
      string curr = atm.tp + "-" + nbAtm.tp;
      if (res.find(curr) == res.end()) {
        curr = nbAtm.tp + "-" + atm.tp;
      }
      if (i) {
        res[curr] += 1;
      } else {
        res[curr] -= 1;
      }
    }
  }
  return res;
}



