#include "gbCnf.h"
#include "KNHome.h"

/*
 * Author: 1mingfei 
 * Date:   2019-05-25
 * Purpose: encoding neighbor atoms of the jumping pairs
 * self-explained
 */

/*
 * Sort points lexicographically 
 * this will only work for the on-lattice cases
 * because it has no buffer range when comparing
 */
inline void sortAtomLexi(vector<KNAtom>& atmList) {
  sort(
    atmList.begin(), atmList.end(),
    [](const KNAtom& a, const KNAtom& b) -> bool
      { return (a.pst[X] < b.pst[X]) || 
               (a.pst[X] == b.pst[X] && a.pst[Y] < b.pst[Y]) ||
               (a.pst[X] == b.pst[X] && a.pst[Y] == b.pst[Y] && 
                a.pst[Z] < b.pst[Z]); }
    );
}

/*
 * Sort input vectors by their first two element
 * e.g. 
 * before sorting:       after sorting:
 * 1 2 4 9               0 2 8 9
 * 2 9 7 3         ==>   0 3 3 4
 * 0 3 3 4         ==>   1 2 4 9
 * 0 2 8 9               2 9 7 3
 */
inline bool sortPairs(const vector<int>& lhs, const vector<int>& rhs) {
  assert(lhs.size() >= 2 && rhs.size() >= 2);
  if (lhs[0] == rhs[0]) {
    return lhs[1] < rhs[1];
  } else {
    return lhs[0] < rhs[0];
  }
}
/* output vectors to file */
inline void writeVector(const string& fname, const vector<int>& v) {
  ofstream ofs(fname, std::ofstream::app);
  for (const auto& val : v) {
    ofs << val << " ";
  }
  ofs << endl;
}

void KNHome::KNEncode() {
  gbCnf cnfModifier(*this);
  vector<string> elems = vsparams["elems"];
  int NConfigs = iparams["NConfigs"];
  int NBars = iparams["NBarriers"];
  double RCut = dparams["RCut"];
  //string fpname = "pairs.txt";//
  string fpname = sparams["PairFile"];
  vector<vector<int>> pairs = readPairs(fpname);
  MPI_Barrier(MPI_COMM_WORLD);
  if (me == 0)
    std::cout << "encoding configurations\n";
  //int quotient = NConfigs / nProcs;
  //int remainder = NConfigs % nProcs;
  int quotient = pairs.size() / nProcs;
  int remainder = pairs.size() % nProcs;
  int nCycle = remainder ? (quotient + 1) : quotient;
  for (int j = 0; j < nCycle; ++j) {
    for (int i = (j * nProcs); i < ((j + 1) * nProcs); ++i) {
      if ((i % nProcs != me) || (i >= NConfigs)) continue;
      //string fname = "in.cfg";//
      string fname = "config" + to_string(pairs[i][0]) + "/s/start.cfg";

      Config cfg = cnfModifier.readCfg(fname);
      vector<int> pair = {pairs[i][2], pairs[i][3]};
      vector<int> tmp = cnfModifier.encodeConfig(cfg, pair, RCut);
      writeVector(to_string(i) + ".txt", tmp);
#ifdef DEBUG
      Config cNew = cfg;
      cNew.atoms.clear();
      cNew.natoms = tmp.size();
      for (unsigned int i = 0 ; i < tmp.size(); ++i) {
        cNew.atoms.push_back(cfg.atoms[tmp[i]]);
      }
      cnfModifier.writeCfgData(cNew, to_string(i) + "_encode.cfg");
#endif
    }
  }

}


vector<vector<int>> KNHome::readPairs(const string& fname) {
  ifstream ifs(fname.empty() ? sparams["datafile"] : fname, std::ifstream::in);
  vector<vector<int>> res;
  vector<string> s;
  string buff;
  getline(ifs, buff); //for comment line Config Barrier atom1 atom2
  while (getline(ifs, buff)) {
    s.clear();
    split(buff, " ", s);
    if (s[0] == "Delta") //for the very last line
      continue;
    vector<int> tmp;
    //assert(s.size() == 4);
    tmp.push_back(stoi(s[1]));
    tmp.push_back(stoi(s[3]));
    for (int i = 5; i < 7; ++i) {
      tmp.push_back(stoi(s[i]));
    }
    res.push_back(tmp);
  }
#ifdef DEBUG
  for (int i = 0; i < res.size(); ++i) {
    writeVector("before.txt", res[i]);
  }
#endif
  sort(res.begin(), res.end(), sortPairs);
#ifdef DEBUG
  for (int i = 0; i < res.size(); ++i) {
    writeVector("after.txt", res[i]);
  }
#endif
  return res;
}

vector<int> KNHome::gbCnf::encodeConfig(Config& cnf,
                                        const vector<int> pair, \
                                        double RCut) {
  assert(pair.size() == 2); //the size of input pair must equals 2
  getNBL(cnf, RCut);
  /* reverse wrap back so the lexi order is always not affected by Periodic 
   * boundary conditions */
  vector<KNAtom> atmList;
  vector<int> res;
  set<int> hsh;
  for (int i = 0; i < 2; ++i) {
    KNAtom atm = cnf.atoms[pair[i]];
    for (const int ii : atm.NBL) {
      KNAtom nbAtm = cnf.atoms[ii];
      for (int j = 0; j < DIM; ++j) {
        if (std::abs(atm.prl[j] - 1.0) < RCut/cnf.length[j]) {
          if ((1.0 - nbAtm.prl[j]) < RCut/cnf.length[j]) {
            nbAtm.prl[j] -= 1.0;
          } else if ((nbAtm.prl[j] - 1.0) < RCut/cnf.length[j]) {
            nbAtm.prl[j] += 1.0;
          }
        }
      }
      if (hsh.find(nbAtm.id) == hsh.end()) {
        hsh.insert(nbAtm.id);
        if ((nbAtm.id != cnf.atoms[pair[0]].id) && \
            (nbAtm.id != cnf.atoms[pair[1]].id)) {
          atmList.push_back(nbAtm);
        }
      }
    }
  }
  sortAtomLexi(atmList);
  for (const auto& atm : atmList) {
    res.push_back(atm.id);
  }
  return res;
}
