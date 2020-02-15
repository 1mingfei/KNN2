//
// Created by Zhucong Xi on 1/28/20.
//
#include "gbCnf.h"
#include "KNHome.h"

/*
Config gbCnf::removeAtoms(Config &inCnf,
                          const string &atomType,
                          unordered_set<int> &atomID) {
  getNBL(inCnf, 3.5);
  Config cnf;
  cnf.cell = inCnf.cell;
  cnf.length = inCnf.length;
  cnf.bvx = inCnf.bvx;
  cnf.tvx = inCnf.tvx;
  cnf.bvy = inCnf.bvy;
  cnf.tvy = inCnf.tvy;
  cnf.bvz = inCnf.bvz;
  cnf.tvz = inCnf.tvz;
  cnf.vacList = inCnf.vacList;
  cnf.ntypes = inCnf.ntypes - 1;
  for (auto &atm : inCnf.atoms) {
    if (atm.tp != atomType) {
      atomID.insert(atm.id);
      cnf.atoms.push_back(atm);
      cnf.natoms++;
    }
  }
  return cnf;
}
*/

unordered_set<int> gbCnf::findSoluteAtoms(const Config &inCnf,
                                          const string &solventAtomType) {
  unordered_set<int> soluteAtomID;
  for (const auto &atm : inCnf.atoms) {
    if (atm.tp != solventAtomType) {
      soluteAtomID.insert(atm.id);
    }
  }
  return soluteAtomID;
}

int gbCnf::helperBFS(const Config &inCnf,
                     const unordered_set<int> &soluteAtomID,
                     unordered_multimap<int, int> &clt2Atm,
                     map<int, int> &atm2Clt) {
//  bool *visited = new bool[soluteAtomID.size()];
//  for (int i = 0; i < soluteAtomID.size(); i++) { visited[i] = false; }

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

      for (int fnnID:inCnf.atoms[atmID].FNNL) {
        // if inFNN in the unvisited set
        it = unvisited.find(fnnID);
        if (it != unvisited.end()) {
          unvisited.erase(it);
          visitQueue.push(fnnID);
        }
      }
    }
    cltID++;
  }
  return cltID;
}

void gbCnf::getLargestClts(const int &numClustersFound,
                           const int &numClustersKept,
                           unordered_multimap<int, int> &clt2Atm,
                           map<int, int> &atm2Clt) {
  vector<vector<int>> keyValueNumMat;
  keyValueNumMat.resize(numClustersFound);
  for (int i = 0; i < numClustersFound; i++) {
    keyValueNumMat[i].push_back(i);
    keyValueNumMat[i].push_back(clt2Atm.count(i));
  }
  //compare second col
  int col = 1;
  sort(keyValueNumMat.begin(),
       keyValueNumMat.end(),
       [col](const vector<int> &lhs, const vector<int> &rhs) -> bool {
         return lhs[col] > rhs[col];
       });

  unordered_multimap<int, int> clt2Atm2;
  map<int, int> atm2Clt2;
  for (int i = 0; i < numClustersKept; i++) {
    int key = keyValueNumMat[i][0];
    auto beg = clt2Atm.equal_range(key).first;
    auto end = clt2Atm.equal_range(key).second;
    for (auto m = beg; m != end; m++) {
      clt2Atm2.insert(pair<int, int>(m->first, m->second));
      atm2Clt2.insert(pair<int, int>(m->second, m->first));
    }
  }
  clt2Atm.clear();
  clt2Atm = clt2Atm2;
  atm2Clt.clear();
  atm2Clt = atm2Clt2;
}
// add FNNs back
void gbCnf::helperAddFNNs(const Config &cnfReference,
                          unordered_multimap<int, int> &clt2Atm,
                          map<int, int> &atm2Clt) {
  for (pair<int, int> i : atm2Clt) {
    int atom = i.first;
    int cluster = i.second;
    for (int j : cnfReference.atoms[atom].FNNL) {
      atm2Clt.insert(pair<int, int>(j, cluster));
    }
  }

  clt2Atm.clear();
  for (pair<int, int> i : atm2Clt) {
    int atom = i.first;
    int cluster = i.second;
    clt2Atm.insert(pair<int, int>(cluster, atom));
  }
}

map<int, int> gbCnf::findAtm2Clts(Config &inCnf,
                                  const int &numClustersKept,
                                  const string &solventAtomType) {
  MPI_Barrier(MPI_COMM_WORLD);
  getNBL(inCnf, 3.5);

  if (me == 0) {
    unordered_set<int> soluteAtomID = findSoluteAtoms(inCnf, solventAtomType);
    unordered_multimap<int, int> clt2Atm;
    map<int, int> atm2Clt;
    int numClustersFound = helperBFS(inCnf, soluteAtomID, clt2Atm, atm2Clt);
    getLargestClts(numClustersFound, numClustersKept, clt2Atm, atm2Clt);
    helperAddFNNs(inCnf, clt2Atm, atm2Clt);
    return atm2Clt;
  } else {
    return map<int, int> {};
  }
}

void KNHome::findClts(gbCnf &inGbCnf) {
  Config inCnf = inGbCnf.readCfg(sparams["initconfig"]);
  map<int, int> atm2Clt = inGbCnf.findAtm2Clts(inCnf,
                                               iparams["numClustersKept"],
                                               sparams["solventAtomType"]);
  if (me == 0) {
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

    vector<int> tmp;
    for (const auto &iter : atm2Clt) {
      outCnf.atoms.push_back(inCnf.atoms[iter.first]);
      tmp.push_back(iter.second);
    }

    inGbCnf.writeCfgAux(outCnf, tmp, "cluster_with_id.cfg");
  }

}