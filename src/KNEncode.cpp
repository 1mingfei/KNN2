#include "gbCnf.h"
#include "KNHome.h"

/* Sort points lexicographically 
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

void KNHome::KNEncode() {
  gbCnf cnfModifier(*this);
  vector<string> elems = vsparams["elems"];
  int NConfigs = iparams["NConfigs"];
  int NBars = iparams["NBarriers"];
  double RCut = dparams["RCut"];
  string fpname = "pairs.txt";//
  vector<pair<int, int>> pairs = readPairs(fpname);
  MPI_Barrier(MPI_COMM_WORLD);
  if (me == 0)
    std::cout << "encoding configurations\n";
  int quotient = NConfigs / nProcs;
  int remainder = NConfigs % nProcs;
  int nCycle = remainder ? (quotient + 1) : quotient;
  for (int j = 0; j < nCycle; ++j) {
    for (int i = (j * nProcs); i < ((j + 1) * nProcs); ++i) {
      if ((i % nProcs != me) || (i >= NConfigs)) continue;
      string fname = ""; //
      cnfModifier.readCfg(fname);
      //pair =  //need sort incoming pairs
      cnfModifier.encodeConfig(cnf, pair, RCut);

    }

}

vector<pair<int, int>> KNHome::readPairs(const string& fname) {
  ifstream ifs(fname.empty() ? sparams["datafile"] : fname, std::ifstream::in);
  vector<pair<int, int>> res;
  vector<vector<int>> inputs;

  vector<string> s;
  string buff;
  getline(ifs, buff);
  sscanf(buff.c_str(), "# %lf %lf", &cnf.oldEngy, &cnf.engy);

  getline(ifs, buff, ' ');
  cnf.natoms = stoi(buff);
  getline(ifs, buff);
  while (getline(ifs, buff)) {
    s.clear();
    split(buff, " ", s);
  }



}

vector<int> KNHome::gbCnf::encodeConfig(const Config& cnf,
                                        const vector<int> pair, \
                                        double RCut) {
  assert(pair.size() == 2); //the size of input pair must equals 2
  getNBL(cnf, Rcut);
  /* reverse wrap back so the lexi order is always not affected by Periodic 
   * boundary conditions */
  vector<KNAtom> atmList;
  vector<int> res;
  set<int> hsh;
  for (int i = 0; i < 2; ++i) {
    KNAtom atm = cnf.atoms[pair[i]];
    for (auto&& nbAtm : atm.NBL) {
      for (int j = 0; j < DIM; ++j) {
        if (std::abs(atm.pst[j] - cnf.length[j]) < Rcut) {
          if ((cnf.length[j] - nbAtm.pst[j]) < Rcut) {
            nbAtm.pst[j] -= cnf.length[j];
          } else if ((nbAtm.pst[j] - cnf.length[j]) < Rcut) {
            nbAtm.pst[j] += cnf.length[j];
          }
        }
      }
      if (hsh.find(nbAtm.id) == hsh.end()) {
        hsh.insert(nbAtm.id);
        if (nbAtm.id != cnf.atoms[0].id) {
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
