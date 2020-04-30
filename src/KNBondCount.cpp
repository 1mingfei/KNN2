#include "gbCnf.h"
#include "KNHome.h"

/* output hashmap to file */
template<class T, class Y>
inline void writeHashmap(const string& fname, \
                         const int i, \
                         const int j, \
                         const map<T, Y>& v, \
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

template<class T, class Y>
inline void writeUnordered_map(const string& fname, \
                               const int i, \
                               const int j, \
                               unordered_map<T, Y>& v, \
                               vector<string>& seq) {
  ofstream ofs(fname, std::ofstream::app);
  ofs << "config " << i << " end " << j << " ";
  for (auto& val : seq)
    ofs << std::setw(4) << v[val] << " ";
  ofs << endl;
}

template<class T, class Y>
inline void writeUnordered_map(const string& fname, \
                               const string& fn, \
                               unordered_map<T, Y>& v, \
                               vector<string>& seq) {
  ofstream ofs(fname, std::ofstream::app);
  ofs << fn << " ";
  for (auto& val : seq)
    ofs << std::setw(4) << v[val] << " ";
  ofs << endl;
}

void KNHome::KNBondCount(gbCnf& cnfModifier) {
  vector<string> elems = vsparams["elems"];
  double RCut = dparams["RCut"];
  string fpname = sparams["PairFile"];
  vector<vector<int>> pairs = readPairs(fpname);

  map<string, int> res;
  res["Al-Al"] = 0;
  res["Al-Mg"] = 0;
  res["Al-Zn"] = 0;
  res["Mg-Mg"] = 0;
  res["Mg-Zn"] = 0;
  res["Zn-Zn"] = 0;
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
      string fname = "config" + to_string(pairs[i][0]) + "/s/start.cfg";

      Config cfg = cnfModifier.readCfg(fname);
      vector<int> pair = {pairs[i][2], pairs[i][3]};
      cnfModifier.getNBL(cfg, RCut);
      map<string, int> pairCount;
      pairCount = cnfModifier.countPairs(cfg, elems, pair);
      writeHashmap("bondCount.txt", pairs[i][0], pairs[i][1], pairCount, \
                   false);
    }
  }
}

map<string, int> gbCnf::countPairs(Config& cnf, \
                                   const vector<string>& elems, \
                                   const vector<int>& pair) {
  map<string, int> res;
  res["Al-Al"] = 0;
  res["Al-Mg"] = 0;
  res["Al-Zn"] = 0;
  res["Mg-Mg"] = 0;
  res["Mg-Zn"] = 0;
  res["Zn-Zn"] = 0;

  if (cnf.atoms[pair[0]].tp == "X" && cnf.atoms[pair[1]].tp == "X")
    return res;


  KNAtom initAtm;
  KNAtom finalAtm;
  if (cnf.atoms[pair[0]].tp != "X") {
    initAtm = cnf.atoms[pair[0]]; //not Vacancy
    finalAtm = cnf.atoms[pair[1]]; // Vacancy
  } else {
    initAtm = cnf.atoms[pair[1]];
    finalAtm = cnf.atoms[pair[0]];
  }

  cout << "init " << initAtm.tp << endl;

  for (const int ii : initAtm.NBL) {
    KNAtom nbAtm = cnf.atoms[ii];
    cout << ii << endl;
    if (nbAtm.tp != "X") {
      string curr = initAtm.tp + "-" + nbAtm.tp;
      if (res.find(curr) == res.end())
        curr = nbAtm.tp + "-" + initAtm.tp;
      res[curr] -= 1;
    }
  }
  cout << " final " << finalAtm.tp << endl;

  for (const int ii : finalAtm.NBL) {
    cout << ii << endl;
    KNAtom nbAtm = cnf.atoms[ii];
    if (nbAtm.tp != "X") {
      string curr = initAtm.tp + "-" + nbAtm.tp;
      if (res.find(curr) == res.end())
        curr = nbAtm.tp + "-" + initAtm.tp;
      res[curr] += 1;
    }
  }

  return res;
}

void KNHome::KNBondCountAll(gbCnf& cnfModifier) {
  vector<string> elems = vsparams["elems"];
  double RCut = dparams["RCut"];
  string fpname = sparams["PairFile"];
  vector<vector<int>> pairs = readPairs(fpname);

  vector<string> outputSeq = {"Al-Al", "Al-Mg", "Al-Zn", \
                              "Mg-Mg", "Mg-Zn", "Zn-Zn"};

  ofstream ofsBC("bondCountAll.txt", std::ofstream::app);

  MPI_Barrier(MPI_COMM_WORLD);
  if (me == 0) {
    std::cout << "count bondings of the configurations\n";
    for (const auto& val : outputSeq)
      ofsBC << std::setw(10) << val << " ";
    ofsBC << "\n";
  }

  int quotient = pairs.size() / nProcs;
  int remainder = pairs.size() % nProcs;
  int nCycle = remainder ? (quotient + 1) : quotient;

  for (int j = 0; j < nCycle; ++j) {
    for (int i = (j * nProcs); i < ((j + 1) * nProcs); ++i) {
      if ((i % nProcs != me) || (i >= pairs.size())) continue;
      string fname = "config" + to_string(pairs[i][0]) \
                     + "/f_" + to_string(pairs[i][1]) + "/end.cfg";

      Config cfg = cnfModifier.readCfg(fname);

      cnfModifier.getNBL(cfg, RCut);
      unordered_map<string, int> pairCount;
      pairCount = cnfModifier.bondCountAll(cfg);
      writeUnordered_map("bondCountAll.txt", \
                         pairs[i][0], pairs[i][1], pairCount, outputSeq);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

void KNHome::KNBondCountList(gbCnf& cnfModifier) {
  vector<string> elems = vsparams["elems"];

  double RCut = dparams["RCut"];
  string fpname = sparams["ListFile"];
  string format = sparams["format"];
  vector<string> outputSeq = {"Al-Al", "Al-Mg", "Al-Zn", \
                              "Mg-Mg", "Mg-Zn", "Zn-Zn"};

  ofstream ofsBC("bondCountAll.txt", std::ofstream::app);
  MPI_Barrier(MPI_COMM_WORLD);
  if (me == 0) {
    std::cout << "count bondings of the configurations\n";
    for (const auto& val : outputSeq)
      ofsBC << std::setw(10) << val << " ";
    ofsBC << "\n";
    vector<string> FileList;
    ifstream file(fpname);
    string str;
    while (getline(file, str)) {
      FileList.push_back(str);
    }
    for (const auto fn : FileList) {
      Config cfg;
      if (format == "vasp")
        cfg = std::move(cnfModifier.readPOSCAR(fn));
      else if (format == "cfg")
        cfg = std::move(cnfModifier.readCfg(fn));
      cnfModifier.wrapAtomPrl(cfg);
      cnfModifier.getNBL_serial(cfg, RCut);
      unordered_map<string, int> pairCount = cnfModifier.bondCountAll(cfg);
      writeUnordered_map("bondCountAll.txt", fn, pairCount, outputSeq);
      pairCount.clear();
    }
  }
}

unordered_map<string, int> gbCnf::bondCountAll(const Config& cfg) {
  unordered_map<string, int> bondsCount;
  string tp1, tp2, bond;
  // int count = 0;
  for (const auto& atom : cfg.atoms) {
    tp1 = atom.tp;
    if (tp1 == "X")
      continue;
    for (const auto& atom2ID : atom.NBL) {
      // if (atom2ID == -1)
      //   continue;
      tp2 = cfg.atoms[atom2ID].tp;
      // cout << (++count) << " " << tp1 << " " << tp2 << " ";
      if (tp2 == "X")
        continue;
      if (tp1.compare(tp2) < 0) {
        bond = tp1;
        bond += "-";
        bond += tp2;
      } else {
        bond = tp2;
        bond += "-";
        bond += tp1;
      }
      // cout << atom.id << " " << atom2ID << " " << bond << " ";

      bondsCount[bond]++;

      // for (auto i : bondsCount)
      //   cout << i.first << " " << i.second << " ";
      // cout << "\n";

    }
  }
  for (auto &[bond, count] : bondsCount) {
    count /= 2;
  }


  return bondsCount;
}