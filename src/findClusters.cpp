#include "gbCnf.h"
#include "KNHome.h"

#define Z 12.0
inline vector<unordered_set<int>> filter( \
                                  const vector<unordered_set<int>>& inMap, \
                                  const int& low, \
                                  const int& high) {
  vector<unordered_set<int>> res = inMap;
  for (int i = 0; i < inMap.size(); ++i) {
    if (low > inMap[i].size() || inMap[i].size() > high)
      res[i].clear();
  }
  return res;
}

inline void sortAtomLexi(vector<KNAtom>& atmList) {
  sort(
    atmList.begin(), atmList.end(),
    [](const KNAtom& a, const KNAtom& b) -> bool
      { return (a.prl[0] < b.prl[0]) ||
               (a.prl[0] == b.prl[0] && a.prl[1] < b.prl[1]) ||
               (a.prl[0] == b.prl[0] && a.prl[1] == b.prl[1] &&
                a.prl[2] < b.prl[2]); }
    );
}

inline bool sizeMatch(const Config& inCnf, \
                      const unordered_set<int>& s, \
                      const Config& refCnf, \
                      vector<double>& lowerLimits, \
                      const double& LC) {
  double Xlimit = refCnf.bvx[0];
  double Ylimit = refCnf.bvy[1];
  double Zlimit = refCnf.bvz[2];
  double xmin = 2.0 * Xlimit, ymin = 2.0 * Ylimit, zmin = 2.0 * Zlimit;
  double xmax = -10.0, ymax = -10.0, zmax = -10.0;
  vector<KNAtom> atmList;
  for (const auto& i : s) {
    KNAtom atm = inCnf.atoms[i];
    xmin = std::min(atm.pst[0], xmin);
    ymin = std::min(atm.pst[1], ymin);
    zmin = std::min(atm.pst[2], zmin);
    xmax = std::max(atm.pst[0], xmax);
    ymax = std::max(atm.pst[1], ymax);
    zmax = std::max(atm.pst[2], zmax);
    atmList.push_back(atm);
  }
  sortAtomLexi(atmList);
  // lowerLimits = {xmin, ymin, zmin};
  lowerLimits = {atmList[0].pst[0], atmList[0].pst[1], atmList[0].pst[2]};

  if (xmax - xmin - 0.5 * LC >= Xlimit)
    return false;
  if (ymax - ymin - 0.5 * LC >= Ylimit)
    return false;
  if (zmax - zmin - 0.5 * LC >= Zlimit)
    return false;
  return true;
}

pair<unordered_set<int>, unordered_set<int>> gbCnf::findSoluteAtoms( \
                                          const Config& inCnf, \
                                          const string& solventAtomType) {
  unordered_set<int> soluteAtomID;
  unordered_set<int> solventAtomID;
  for (const auto& atm : inCnf.atoms) {
    if (atm.tp == solventAtomType || atm.tp == "Xe") { // Xe always like this
      solventAtomID.insert(atm.id);
    } else {
      soluteAtomID.insert(atm.id);
    }
  }
  return make_pair(soluteAtomID, solventAtomID);
}

int gbCnf::helperBFS(const Config& inCnf, \
                     const unordered_set<int>& soluteAtomID, \
                     unordered_multimap<int, int>& clt2Atm, \
                     map<int, int>& atm2Clt) {
  unordered_set<int> unvisited = soluteAtomID;
  queue<int> visitQueue;

  int cltID = 0;
  while (!unvisited.empty()) {
    // Find the next element
    auto it = unvisited.begin();
    int atmID = *it;
    visitQueue.push(atmID);
    unvisited.erase(it);

    while (!visitQueue.empty()) {
      atmID = visitQueue.front();
      visitQueue.pop();
      // Add to maps
      clt2Atm.insert(pair<int, int>(cltID, atmID));
      atm2Clt.insert(pair<int, int>(atmID, cltID));

      for (int fnnID : inCnf.atoms[atmID].FNNL) {
        if (fnnID == -1)
          continue;
        // if inFNN in the unvisited set
        it = unvisited.find(fnnID);
        if (it != unvisited.end()) {
          unvisited.erase(it);
          visitQueue.push(fnnID);
        }
      }
    }
    ++cltID;
  }
  return cltID;
}

int gbCnf::helperBFSRmMtrx(const Config& inCnf, \
                           unordered_multimap<int, int>& clt2Atm, \
                           map<int, int>& atm2Clt, \
                           const int& index, \
                           const string& solventAtomType, \
                           int& numAtomsLeft) {

  unordered_set<int> visited;
  queue<int> Q;

  Q.push(index);
  // search for connected solvent elements
  while (!Q.empty()) {
    int size = Q.size();
    for (int i = 0; i < size; ++i) {

      auto atmID = Q.front();
      Q.pop();
      if (visited.find(atmID) != visited.end())
        continue;
      visited.insert(atmID);

      for (int fnnID : inCnf.atoms[atmID].FNNL) {
        if (fnnID == -1)
          continue;
        if ((inCnf.atoms[fnnID].tp == "Xe") \
            || (inCnf.atoms[fnnID].tp == solventAtomType)) {
          Q.push(fnnID);
        }
      }
    }
  }

  numAtomsLeft -= visited.size();

  int cltID = 0;
  unordered_set<int> unvisited;
  queue<int> visitQueue;

  for (const auto& atm : inCnf.atoms) {
    if (visited.find(atm.id) == visited.end())
      unvisited.insert(atm.id);
  }

  while (!unvisited.empty()) {
    // Find the next element
    auto it = unvisited.begin();
    int atmID = *it;
    visitQueue.push(atmID);
    unvisited.erase(it);

    while (!visitQueue.empty()) {
      atmID = visitQueue.front();
      visitQueue.pop();
      // Add to maps
      clt2Atm.insert(pair<int, int>(cltID, atmID));
      atm2Clt.insert(pair<int, int>(atmID, cltID));

      for (int fnnID : inCnf.atoms[atmID].FNNL) {
        if (fnnID == -1)
          continue;
        // if inFNN in the unvisited set
        it = unvisited.find(fnnID);
        if (it != unvisited.end()) {
          unvisited.erase(it);
          visitQueue.push(fnnID);
        }
      }
    }
    ++cltID;
  }
  return cltID;
}

int gbCnf::getLargestClts(const int& numClustersFound, \
                          const int& numClustersKept, \
                          unordered_multimap<int, int>& clt2Atm, \
                          map<int, int>& atm2Clt, \
                          const int& smallest_cluster) {

  vector<vector<int>> keyValueNumMat;
  keyValueNumMat.resize(numClustersFound);
  for (int i = 0; i < numClustersFound; ++i) {
    keyValueNumMat[i].push_back(i);
    keyValueNumMat[i].push_back(clt2Atm.count(i));
  }
  //compare second col
  int col = 1;
  sort(keyValueNumMat.begin(),
       keyValueNumMat.end(),
       [col](const vector<int>& lhs, const vector<int>& rhs) -> bool {
         return lhs[col] > rhs[col];
       });

  unordered_multimap<int, int> clt2Atm2;
  map<int, int> atm2Clt2;
  int safeNumCluster = numClustersKept < numClustersFound \
                        ? numClustersKept : numClustersFound;
  for (int i = 0; i < safeNumCluster; ++i) {
    int key = keyValueNumMat[i][0];
    auto beg = clt2Atm.equal_range(key).first;
    auto end = clt2Atm.equal_range(key).second;
    int dist = distance(beg, end);
    if (dist > smallest_cluster) {
      for (auto m = beg; m != end; ++m) {
        clt2Atm2.insert(pair<int, int>(i, m->second));
        atm2Clt2.insert(pair<int, int>(m->second, i));
      }
    }
  }
  clt2Atm.clear();
  clt2Atm = clt2Atm2;
  atm2Clt.clear();
  atm2Clt = atm2Clt2;

  return safeNumCluster;
}

bool gbCnf::validSolventCluster(const Config& cnf, \
                                const int& i, \
                                const string& solventAtomType, \
                                const int& solventBoudCriteria, \
                                const unordered_set<int>& atm_other) {
  int count = 0;
  for (int j : cnf.atoms[i].FNNL) {
    if (j == -1)
      continue;
    if ((cnf.atoms[j].tp != solventAtomType) \
        && (cnf.atoms[j].tp != "Xe") \
        && (cnf.atoms[j].tp != "X") \
        && (atm_other.find(j) != atm_other.end()))
      ++count;
  }
  return (count >= solventBoudCriteria);
}

// add FNNs back
void gbCnf::helperAddFNNs(const Config& cnfReference, \
                          unordered_multimap<int, int>& clt2Atm, \
                          map<int, int>& atm2Clt, \
                          const string& solventAtomType, \
                          const int& solventBoudCriteria, \
                          const int& numClustersFound) {
  unordered_set<int> atm_other;
  for (pair<int, int> i : atm2Clt) {
    atm_other.insert(i.first);
  }
  for (pair<int, int> i : atm2Clt) {
    int atom = i.first;
    int cluster = i.second;
    for (int j : cnfReference.atoms[atom].FNNL) {
      if (j == -1)
        continue;
      if (validSolventCluster(cnfReference, j, solventAtomType, \
                              solventBoudCriteria, atm_other)) {
        atm2Clt.insert(pair<int, int>(j, cluster));
      }
    }
  }

  clt2Atm.clear();
  for (pair<int, int> i : atm2Clt) {
    int atom = i.first;
    int cluster = i.second;
    clt2Atm.insert(pair<int, int>(cluster, atom));
  }

#ifdef OUTPUT_CLT_DISTRI
  vector<vector<int>> keyValueNumMat;
  keyValueNumMat.resize(numClustersFound);
  for (int i = 0; i < numClustersFound; ++i) {
    keyValueNumMat[i].push_back(i);
    keyValueNumMat[i].push_back(clt2Atm.count(i));
  }

  int col = 1;
  sort(keyValueNumMat.begin(),
       keyValueNumMat.end(),
       [col](const vector<int>& lhs, const vector<int>& rhs) -> bool {
         return lhs[col] > rhs[col];
       });


  map<int, int> mp;
  for (int i = 0; i < numClustersFound; ++i) {
    int key = keyValueNumMat[i][0];
    auto beg = clt2Atm.equal_range(key).first;
    auto end = clt2Atm.equal_range(key).second;
    int dist = distance(beg, end);
    ++mp[dist];
  }
  for (const auto& it : mp) {
    cout << it.first << " " << it.second << "\n";
  }
  cout << "\n";
#endif
}

map<int, int> gbCnf::findAtm2Clts(Config& inCnf, \
                                  const int& numClustersKept, \
                                  const string& solventAtomType, \
                                  const int& solventBoudCriteria, \
                                  const int& smallest_cluster) {
  MPI_Barrier(MPI_COMM_WORLD);
  if (nProcs == 1)
    getNBL_serial(inCnf, 3.5);
  else
    getNBL(inCnf, 3.5);

  if (me == 0) {
    auto soluteAtomID = findSoluteAtoms(inCnf, solventAtomType).first;
    unordered_multimap<int, int> clt2Atm;
    map<int, int> atm2Clt;
    int numClustersFound = helperBFS(inCnf, soluteAtomID, clt2Atm, atm2Clt);
    int safe = getLargestClts(numClustersFound, numClustersKept, clt2Atm, \
                              atm2Clt, smallest_cluster);
    helperAddFNNs(inCnf, clt2Atm, atm2Clt, \
                  solventAtomType, solventBoudCriteria, safe);
    return atm2Clt;
  } else {
    return map<int, int>{};
  }
}

map<int, int> gbCnf::findAtm2CltsRmMtrx(Config& inCnf, \
                                        const string& solventAtomType, \
                                        int& numAtomsLeft, \
                                        const int& smallest_cluster) {
  MPI_Barrier(MPI_COMM_WORLD);

  if (nProcs == 1)
    getNBL_serial(inCnf, 3.5);
  else
    getNBL(inCnf, 3.5);

  auto ssID = findSoluteAtoms(inCnf, solventAtomType);
  auto soluteAtomID = ssID.first;
  auto solventAtomID = ssID.second;

  vector<int> solventAtomID_v;
  for (const auto& elem : solventAtomID) {
    solventAtomID_v.push_back(elem);
  }
  unordered_multimap<int, int> clt2Atm;
  map<int, int> atm2Clt;

  int index = (rand() + me) % solventAtomID_v.size();

  int numClustersFound = helperBFSRmMtrx(inCnf, \
                                         clt2Atm, \
                                         atm2Clt, \
                                         solventAtomID_v[index], \
                                         solventAtomType, \
                                         numAtomsLeft);

  // cout << "# " << me << " " << solventAtomID_v[index] << " " \
  //      << numClustersFound << " " << numAtomsLeft << " " \
  //      << atm2Clt.size() << endl;

  int tmp = getLargestClts(numClustersFound, 108000, \
                           clt2Atm, atm2Clt, smallest_cluster);

  return atm2Clt;
}

void KNHome::findClts(gbCnf& inGbCnf, \
                      const string& fname, \
                      const string& mode) {
  Config inCnf = inGbCnf.readCfg(fname);
  map<int, int> atm2Clt;

  int numAtomsLeft = inCnf.natoms;

  int smallest_cluster = (iparams["smallestClusterCriteria"] == 0) \
                          ? 3 : iparams["smallestClusterCriteria"];

  if (mode == "clusterCount") {
    atm2Clt = inGbCnf.findAtm2Clts(inCnf, \
                                   iparams["numClustersKept"], \
                                   sparams["solventAtomType"],
                                   iparams["solventBoudCriteria"], \
                                   smallest_cluster);
  } else if (mode == "clusterCountRemoveMatrix") {
    atm2Clt = inGbCnf.findAtm2CltsRmMtrx(inCnf, \
                                         sparams["solventAtomType"], \
                                         numAtomsLeft, \
                                         smallest_cluster);
  }

  int pick = 0, maxVal = -1;
  MPI_Barrier(MPI_COMM_WORLD);
  if (nProcs > 1) {
    int* buffData = new int [nProcs];
    MPI_Allgather(&numAtomsLeft, 1, MPI_INT, buffData, 1, MPI_INT, \
                  MPI_COMM_WORLD);

    for (int i = 0; i < nProcs; ++i) {
      if (buffData[i] > maxVal) {
        pick = i;
        maxVal = buffData[i];
      }
    }

    delete[] buffData;
  }

  if (me == pick) {
    Config outCnf;
    outCnf.cell = inCnf.cell;
    outCnf.length = inCnf.length;
    outCnf.bvx = inCnf.bvx;
    outCnf.tvx = inCnf.tvx;
    outCnf.bvy = inCnf.bvy;
    outCnf.tvy = inCnf.tvy;
    outCnf.bvz = inCnf.bvz;
    outCnf.tvz = inCnf.tvz;
    outCnf.vacList = inCnf.vacList;
    outCnf.natoms = atm2Clt.size();

    vector<int> cltId;
    for (const auto& item : atm2Clt) {
      outCnf.atoms.push_back(inCnf.atoms[item.first]);
      cltId.push_back(item.second);
    }

    vector<string> str;
    split(fname, ".", str);
    string oFName;
    oFName = str[0] + "_cluster.cfg";
    inGbCnf.writeCfgAux(outCnf, cltId, oFName);

    ofs.open("clusters_info.txt", std::ofstream::out | std::ofstream::app);
    std::map<string, int> names;
    vector<string> elems = vsparams["elems"];
    for (const auto& elem : elems) {
      names.insert(make_pair(elem, 0));
    }
    for (const auto& atm : outCnf.atoms) {
      ++names[atm.tp];
    }
    ofs << str[0] << " ";
    for (auto& name : names) {
      if (name.first == "X") continue;
      ofs << name.second << " ";
    }
    ofs << endl;
    ofs.close();
  }
}

void KNHome::loopConfigCluster(gbCnf& inGbCnf, const string& mode) {
  long long initNum = (iparams["initNum"] == 0) ? 0 : iparams["initNum"];
  long long increment = (iparams["increment"] == 0) ? 0 : iparams["increment"];
  long long finalNum = (iparams["finalNum"] == 0) ? 0 : iparams["finalNum"];
  for (long long i = initNum; i <= finalNum; i += increment) {
    string fname = to_string(i) + ".cfg";
    calSRO(inGbCnf, fname);
    findClts(inGbCnf, fname, mode);
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

void KNHome::loopConfigSRO(gbCnf& inGbCnf) {
  long long initNum = (iparams["initNum"] == 0) ? 0 : iparams["initNum"];
  long long increment = (iparams["increment"] == 0) ? 0 : iparams["increment"];
  long long finalNum = (iparams["finalNum"] == 0) ? 0 : iparams["finalNum"];
  for (long long i = initNum; i <= finalNum; i += increment) {
    string fname = to_string(i) + ".cfg";
    calSRO(inGbCnf, fname);
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

void KNHome::loopClusterStat(gbCnf& inGbCnf) {
  long long initNum = (iparams["initNum"] == 0) ? 0 : iparams["initNum"];
  long long increment = (iparams["increment"] == 0) ? 0 : iparams["increment"];
  long long finalNum = (iparams["finalNum"] == 0) ? 0 : iparams["finalNum"];
  ofs.open("clusters_size_distr.txt", std::ofstream::app);
  ofs << "filename ";
  for (int i = 0; i < viparams["factors"].size(); ++i) {
    ofs << viparams["factors"][i] << "_count ";
  }
  ofs << ">=" \
      << to_string(viparams["factors"][viparams["factors"].size() - 1]) \
      << "_count ";


  for (int i = 0; i < viparams["factors"].size(); ++i) {
    ofs << viparams["factors"][i] << "_mean " \
        << viparams["factors"][i] << "_std ";
  }
  ofs << ">=" \
      << to_string(viparams["factors"][viparams["factors"].size() - 1]) \
      << "_mean ";
  ofs << ">=" \
      << to_string(viparams["factors"][viparams["factors"].size() - 1]) \
      << "std ";

  for (int i = 0; i < viparams["factors"].size(); ++i) {
    ofs << viparams["factors"][i] << "_Al_conc " \
        << viparams["factors"][i] << "_Mg_conc " \
        << viparams["factors"][i] << "_Zn_conc " ;
  }
  ofs << ">=" \
      << to_string(viparams["factors"][viparams["factors"].size() - 1]) \
      << "_Al_conc ";
  ofs << ">=" \
      << to_string(viparams["factors"][viparams["factors"].size() - 1]) \
      << "_Mg_conc ";
  ofs << ">=" \
      << to_string(viparams["factors"][viparams["factors"].size() - 1]) \
      << "_Zn_conc ";

  for (int i = 0; i < viparams["factors"].size(); ++i) {
    ofs << viparams["factors"][i] << "_Zn/Mg_ratio " ;
    ofs << viparams["factors"][i] << "_Al/(Zn+Mg)_ratio " ;
  }
  ofs << ">=" \
      << to_string(viparams["factors"][viparams["factors"].size() - 1]) \
      << "_Zn/Mg_ratio ";
  ofs << ">=" \
      << to_string(viparams["factors"][viparams["factors"].size() - 1]) \
      << "_Al/(Zn+Mg)_ratio ";

  ofs << "\n";

  for (long long i = initNum; i <= finalNum; i += increment) {
    ofs << i << " ";
    string fname = to_string(i) + "_cluster.cfg";
    clusterStat(inGbCnf, fname);
  }
}

inline void online_variance(const int& x, \
                            map<int, int>& mp, \
                            map<int, pair<double, double>>& mp_mu_std) {

  for (auto&& j : mp) {
    if (x < j.first) {
      ++j.second;
      double delta = static_cast<double>(x - mp_mu_std[j.first].first);
      mp_mu_std[j.first].first += delta / static_cast<double>(j.second);
      double var = delta * static_cast<double>(x - mp_mu_std[j.first].first);
      if (j.second >= 2)
        mp_mu_std[j.first].second = var / static_cast<double>(j.second - 1);
      else
        mp_mu_std[j.first].second = 0.0;
      break;
    }
  }

}

void KNHome::clusterStat(gbCnf& inGbCnf, const string& fname) {
  cout << "working on " << fname << "\n";
  vector<unordered_set<int>> oldmap;
  Config inCnf = inGbCnf.readCfgCluster(fname, oldmap);

  vector<int> interval = viparams["factors"]; // size interval
  map<int, int> mp;
  map<int, pair<double, double>> mp_mu_std; // mean size of cluster & variance
  map<int, vector<int>> mp_spe_count; // species count
  vector<int> nums = viparams["nums"]; // input total number of atoms for spe
  for (int i = 0; i < interval.size(); ++i) {
    mp[interval[i]] = 0;
    mp_mu_std[interval[i]] = {0.0, 0.0};
    mp_spe_count[interval[i]] = vector<int>(3, 0); //Al Mg Zn, Xe treated as Al now
  }

  mp[10000] = 0; // for the last bin
  mp_mu_std[10000] = {0.0, 0.0};
  mp_spe_count[10000] = vector<int>(3, 0);
  for (int i = 0; i < oldmap.size(); ++i) {
    online_variance(oldmap[i].size(), mp, mp_mu_std);
    // count species
    for (auto&& j : mp) {
      if (oldmap[i].size() < j.first) {
        for (auto atomID : oldmap[i]) {
          string type = inCnf.atoms[atomID].tp;
          if (type == "Al" || "Xe")
            ++mp_spe_count[j.first][0];
          if (type == "Mg")
            ++mp_spe_count[j.first][1];
          if (type == "Zn")
            ++mp_spe_count[j.first][2];
        }
        break;
      }
    }
  }

  for (auto i : mp) {
    ofs << i.second << " ";
    cout << i.second << " ";
  }
  for (auto i : mp_mu_std) {
    ofs << i.second.first << " " << std::sqrt(i.second.second) << " ";
  }
  for (auto i : mp_spe_count) {
    for (int j = 0; j < i.second.size(); ++j) {
      ofs << static_cast<double>(i.second[j]) / static_cast<double>(nums[j]) \
          << " ";
    }
  }
  for (auto i : mp_spe_count) {
    ofs << static_cast<double>(i.second[2]) / static_cast<double>(i.second[1])
        << " "
        << static_cast<double>(i.second[0]) / \
           static_cast<double>(i.second[1] + i.second[2]) << " ";

  }
  ofs << "\n";
}

void KNHome::calSRO(gbCnf& inGbCnf, const string& fname) {
  Config inCnf = inGbCnf.readCfg(fname);
  MPI_Barrier(MPI_COMM_WORLD);

  if (nProcs == 1)
    inGbCnf.getNBL_serial(inCnf, 3.5);
  else
    inGbCnf.getNBL(inCnf, 3.5);
  if (me == 0) {
    ofs.open("clusters_SRO.txt", std::ofstream::out | std::ofstream::app);

    vector<string> elems = vsparams["elems"];
    unordered_map<string, int> h_bond;
    unordered_map<string, int> h_elem;

    if (fname == "0.cfg") {
    ofs << "#SRO ";
      for (int i = 0; i < elems.size(); ++i) {
        for (int j = i; j < elems.size(); ++j) {
          if (i != j) {
            ofs << "SRO(" << elems[i] << "-" << elems[j] << ")      ";
            h_bond[elems[i] + "-" + elems[j]] = 0;
          }
        }
      }
      for (int i = 0; i < elems.size(); ++i) {
        for (int j = i; j < elems.size(); ++j) {
          if (i != j) {
            ofs << elems[i] << "-" << elems[j] << "      ";
            h_bond[elems[i] + "-" + elems[j]] = 0;
          }
        }
      }
      ofs << "\n";
    }

    int size = inCnf.natoms;
    double sizeSq = static_cast<double>(size) * static_cast<double>(size);
    for (int i = 0; i < size; ++i) {
      auto&& atm = inCnf.atoms[i];
      ++h_elem[atm.tp];
      for (int fnnID : inCnf.atoms[atm.id].FNNL) {
        if (fnnID == -1)
          continue;
        auto&& fnn = inCnf.atoms[fnnID];
        ++h_bond[atm.tp + "-" + fnn.tp];
      }
    }

    vector<string> str;
    split(fname, ".", str);
    ofs << str[0] << " ";
    for (int i = 0; i < elems.size(); ++i) {
      for (int j = i; j < elems.size(); ++j) {
        if (i != j) {
          double zeta = static_cast<double>(h_bond[elems[i] + "-" + elems[j]])
                        / Z / static_cast<double>(size);
          double res = zeta * static_cast<double>(h_elem[elems[i]]) \
                             * static_cast<double>(h_elem[elems[j]]) / sizeSq;
          ofs << setprecision(9) << (1.0 - res) << " ";
        }
      }
    }
    for (int i = 0; i < elems.size(); ++i) {
      for (int j = i; j < elems.size(); ++j) {
        if (i != j) {
          ofs << h_bond[elems[i] + "-" + elems[j]]  << "  ";
        }
      }
    }
    ofs << "\n";

    ofs.close();
  }
}


void KNHome::resizeCluster(gbCnf& inGbCnf, const string& fname) {
  vector<unordered_set<int>> oldmap;
  Config inCnf = inGbCnf.readCfgCluster(fname, oldmap);
  inGbCnf.cnvprl2pst(inCnf);

  // fliter size
  int low = iparams["low"] ? iparams["low"] : 15;
  int high = iparams["high"] ? iparams["high"] : 100;
  vector<unordered_set<int>> map = filter(oldmap, low, high);
  MPI_Barrier(MPI_COMM_WORLD);

  if (me == 0) {
    vector<int> factors = viparams["factors"];
    double LC = dparams["LC"];
    Config refCnf = inGbCnf.getFCCConv(LC, "Al", factors);
    inGbCnf.cnvprl2pst(refCnf);
    int count = 0;
    int NConfigs = iparams["NConfigs"] < map.size()
                  ? iparams["NConfigs"] : map.size();
    int validCount = 0;
    while(validCount < NConfigs) {
      if (map[validCount++].empty())
        continue;
      vector<double> lowerLimits = {0.0, 0.0, 0.0};
      if (sizeMatch(inCnf, map[validCount], refCnf, lowerLimits, LC)) {
        Config outCnf = inGbCnf.Complete(inCnf, \
                                         map[validCount], \
                                         refCnf, \
                                         lowerLimits);
        string outfname = to_string(count++) + "_" + to_string(validCount) \
                        + ".cfg";
        inGbCnf.writeCfgData(outCnf, outfname);
      }
    }
  }
}

Config gbCnf::Complete(const Config& inCnf, \
                       const unordered_set<int>& s, \
                       const Config& refCnf, \
                       vector<double>& lowerLimits) {
  Config outCnf;
  outCnf.cell = refCnf.cell;
  outCnf.length = refCnf.length;
  outCnf.bvx = refCnf.bvx;
  outCnf.tvx = refCnf.tvx;
  outCnf.bvy = refCnf.bvy;
  outCnf.tvy = refCnf.tvy;
  outCnf.bvz = refCnf.bvz;
  outCnf.tvz = refCnf.tvz;

  for (auto i : s) {
    KNAtom atm = inCnf.atoms[i];
    for (int j = 0; j < 3; ++j) {
      atm.pst[j] -= lowerLimits[j];
      if (atm.pst[j] > outCnf.length[j])
        atm.pst[j] -= outCnf.length[j];
      if (atm.pst[j] < 0.0)
        atm.pst[j] += outCnf.length[j];
    }

    outCnf.atoms.push_back(atm);
  }
  cnvpst2prl(outCnf);
  vector<KNAtom> atomList = outCnf.atoms;
  int count = 0;
  atomList.insert(atomList.end(), refCnf.atoms.begin(), refCnf.atoms.end());


  vector<KNAtom> outList;
  for (int i = 0; i < atomList.size(); ++i) {
    bool isDup = false;
    for (int j = 0; j < outList.size(); ++j) {
      // if (samePos(atomList[i], outList[j]))
      if (calDistPrl(outCnf.length, atomList[i], outList[j]) < 0.2)
        isDup = true;
    }
    if (!isDup)
      outList.push_back(atomList[i]);
  }

  outCnf.atoms = outList;
  outCnf.natoms = outCnf.atoms.size();
  cnvpst2prl(outCnf);
  return outCnf;
}