#include "LSEvent.h"

namespace LS {

LSEvent::LSEvent(const pair<int, int>& inPair)
: KMCEvent(inPair) {}

void LSEvent::exeEvent(Config& cnf, \
                        const double& RCut) {
  int first = jumpPair.first;
  int second = jumpPair.second;

  // first : vac
  // second : other element

  /* 1 */
  // std::swap(cnf.atoms[first].id, cnf.atoms[second].id);
  /* 1 */

  std::swap(cnf.atoms[first].prl, cnf.atoms[second].prl);
  std::swap(cnf.atoms[first].pst, cnf.atoms[second].pst);

  /* 2 jump pair neighbour list
    e.g.
    1  2  3  4  5                    1  2  3  4  5
    6  7  8  9  10          ==>      6 13  8  9 10
    11 12 13 14 15                  11 12  7 14 15
    16 17 18 19 20                  16 17 18 19 20
    left
    NBL:
    7 : 1 2 3 6 8 11 12 13
    13: 8 9 12 14 17 18 19 7
    right
    NBL:
    7 : 8 9 12 14 17 18 19 13
    13: 1 2 3 6 8 11 12 7
  */

  std::swap(cnf.atoms[first].NBL, cnf.atoms[second].NBL);
  std::replace(cnf.atoms[first].NBL.begin(), cnf.atoms[first].NBL.end(), \
          first, second);
  std::replace(cnf.atoms[second].NBL.begin(), cnf.atoms[second].NBL.end(), \
          second, first);

  std::swap(cnf.atoms[first].FNNL, cnf.atoms[second].FNNL);
  std::replace(cnf.atoms[first].FNNL.begin(), cnf.atoms[first].FNNL.end(), \
          first, second);
  std::replace(cnf.atoms[second].FNNL.begin(), cnf.atoms[second].FNNL.end(), \
          second, first);

  /* 3 neighbor atom's neighbor list */
  unordered_set<int> tmpSet;
  for (const int& i : cnf.atoms[first].NBL)
    tmpSet.insert(i);
  for (const int&  i : cnf.atoms[second].NBL)
    tmpSet.insert(i);

  vector<int> tmpList(tmpSet.begin(), tmpSet.end());
  for (int i : tmpList) {
    if (i == first || i == second)
      continue;
    for (int j = 0; j < cnf.atoms[i].NBL.size(); ++j) {
      if (cnf.atoms[i].NBL[j] == first) {
        cnf.atoms[i].NBL[j] = second;
        continue;
      }
      if (cnf.atoms[i].NBL[j] == second)
        cnf.atoms[i].NBL[j] = first;
    }

#ifdef DEBUGJUMP
    cout << "atom " << i << endl;
    for (const auto& nei : cnf.atoms[i].NBL)
      cout << nei << " ";
    cout << endl;
#endif

  }

  tmpSet.clear();
  for (const int& i : cnf.atoms[first].FNNL)
    tmpSet.insert(i);
  for (const int&  i : cnf.atoms[second].FNNL)
    tmpSet.insert(i);
  tmpList.clear();
  for (const int&  i : tmpSet)
    tmpList.emplace_back(i);

  for (int i : tmpList) {
    if (i == first || i == second)
      continue;
    for (int j = 0; j < cnf.atoms[i].FNNL.size(); ++j) {
      if (cnf.atoms[i].FNNL[j] == first) {
        cnf.atoms[i].FNNL[j] = second;
        continue;
      }
      if (cnf.atoms[i].FNNL[j] == second)
        cnf.atoms[i].FNNL[j] = first;
    }

#ifdef DEBUGJUMP
    cout << "atom " << i << endl;
    for (const auto& nei : cnf.atoms[i].FNNL)
      cout << nei << " ";
    cout << endl;
#endif

  }
}

} // end namespace LS