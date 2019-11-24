#include "KMCEvent.h"
// #include "gbCnf.h"

/*calculate distance between one atom in configuration and one from ref*/
inline double calDistPrl(const vector<double> length, \
                         const KNAtom& atm1,\
                         const KNAtom& atm2) {
  double xi = atm1.prl[X];
  double xj = atm2.prl[X];
  double yi = atm1.prl[Y];
  double yj = atm2.prl[Y];
  double zi = atm1.prl[Z];
  double zj = atm2.prl[Z];
  double a, b, c;

  if (xj - xi >= 0.5) {
    a = (xi - xj + 1.0);
  } else if (xj - xi <  -0.5) {
    a = (xi - xj - 1.0); 
  } else {
    a = xi - xj;
  }

  if (yj - yi >= 0.5) {
    b = (yi - yj + 1.0);
  } else if (yj - yi <  -0.5) {
    b = (yi - yj - 1.0); 
  } else {
    b = yi - yj;
  }

  if (zj - zi >= 0.5) {
    c = (zi - zj + 1.0);
  } else if (zj - zi <  -0.5) {
    c = (zi - zj - 1.0); 
  } else {
    c = zi - zj;
  }
  a *= length[X];
  b *= length[Y];
  c *= length[Z];

  double dist = sqrt(a*a + b*b + c*c);
  return dist;
}

KMCEvent::KMCEvent(const pair<int, int>& inPair)
: jumpPair(inPair) {}

KMCEvent::KMCEvent() {
}

KMCEvent::~KMCEvent() {}

double KMCEvent::getRate() const {
  return rate;
}

double KMCEvent::getProb() const {
  return prob;
}

double KMCEvent::getcProb() const {
  return cProb;
}

pair<int, int> KMCEvent::getJumpPair() const {
  return jumpPair;
}


void KMCEvent::calProb(const double& sum) {
  prob = rate / sum;
}

void KMCEvent::setcProb(const double& curr) {
  cProb = curr;
}

void KMCEvent::setRate(const double& inRate) {
  rate = inRate;
}

void KMCEvent::exeEvent(Config& cnf, \
                        unordered_map<int, vector<int>>& jumpList,\
                        const double& RCut) {
  int first = jumpPair.first;
  int second = jumpPair.second;
  /*
    first : vac
    second : other element
    three things to swap here:
    // 1) element id
    2) element coordinates
    3) element neighbour list
    4) neighbor atom's neighbor list
    5) jumpList
   */
  /* 1 */
  // std::swap(cnf.atoms[first].id, cnf.atoms[second].id);
  /* 2 */
  std::swap(cnf.atoms[first].prl, cnf.atoms[second].prl);
  std::swap(cnf.atoms[first].pst, cnf.atoms[second].pst);

  /* 3 
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

#ifdef DEBUGJUMP
    cout << "atom " << first << endl;
    for (const auto& nei : cnf.atoms[first].NBL)
      cout << nei << " ";
    cout << endl;
    cout << "atom " << second << endl;
    for (const auto& nei : cnf.atoms[second].NBL)
      cout << nei << " ";
    cout << endl;
#endif

  std::swap(cnf.atoms[first].NBL, cnf.atoms[second].NBL);
  // replace(cnf.atoms[first].NBL.begin(), cnf.atoms[first].NBL.end(), \
  //         first, second);
  // replace(cnf.atoms[second].NBL.begin(), cnf.atoms[second].NBL.end(), \
  //         second, first);

#ifdef DEBUGJUMP
    cout << "atom " << first << endl;
    for (const auto& nei : cnf.atoms[first].NBL)
      cout << nei << " ";
    cout << endl;
    cout << "atom " << second << endl;
    for (const auto& nei : cnf.atoms[second].NBL)
      cout << nei << " ";
    cout << endl;
#endif

  /* 4 */
  vector<int> tmpList;
  for (int i : cnf.atoms[first].NBL)
    tmpList.push_back(i);
  for (int i : cnf.atoms[second].NBL)
    tmpList.push_back(i);
  for (int i : tmpList) {
    int count = 0;
    if (find(cnf.atoms[i].NBL.begin(), \
             cnf.atoms[i].NBL.end(), \
             first) != cnf.atoms[i].NBL.end()) 
      ++count;
    if (find(cnf.atoms[i].NBL.begin(), \
             cnf.atoms[i].NBL.end(), \
             second) != cnf.atoms[i].NBL.end()) 
      ++count;
    if (count == 2) //nothing need to be done
      continue;
    replace(cnf.atoms[i].NBL.begin(), cnf.atoms[i].NBL.end(), \
            first, -first);
    replace(cnf.atoms[i].NBL.begin(), cnf.atoms[i].NBL.end(), \
            second, first);
    replace(cnf.atoms[i].NBL.begin(), cnf.atoms[i].NBL.end(), \
            -first, second);
    // if (i == first || i == second)
    //   continue;
    // for (int j = 0; j < cnf.atoms[i].NBL.size(); ++j) {
    //   if (cnf.atoms[i].NBL[j] == first)
    //     cnf.atoms[i].NBL[j] = second;
    //   if (cnf.atoms[i].NBL[j] == second)
    //     cnf.atoms[i].NBL[j] = first;
    // }

#ifdef DEBUGJUMP
    cout << "atom " << i << endl;
    for (const auto& nei : cnf.atoms[i].NBL)
      cout << nei << " ";
    cout << endl;
#endif


  }


  /* 5: to do update jumpList by recalculating atoms within cutoff of 3.0A */
  /* find the previous vac; update its neighbor list */
  vector<int> res;
  for (int i = 0; i < cnf.atoms.size(); ++i) {
    if ((i == first) || (cnf.atoms[i].tp == "X"))
      continue;
    double dist = calDistPrl(cnf.length, \
                             cnf.atoms[first], \
                             cnf.atoms[i]);
    if (dist <= RCut)
      res.push_back(i);
  }
  jumpList[first] = std::move(res);
  /* update other jump list by swap vac and the other element */
  for (pair<int, vector<int>> i : jumpList) {
    if (i.first == first)
      continue;
    for (int j = 0; j < i.second.size(); ++j) {
      if (i.second[j] == first)
        i.second[j] = second;
      if (i.second[j] == second)
        i.second[j] = first;
    }
  }

#ifdef DEBUGJUMP
  cout << " new method: " << endl;
  vector<int> vacList;
  for (auto&& it = jumpList.begin(); it != jumpList.end(); ++it)
    vacList.push_back(it->first);
  for (int i = 0; i < vacList.size(); ++i) {
    cout << vacList[i] << " size: " << jumpList[vacList[i]].size() << endl;
    for (int j = 0; j < jumpList[vacList[i]].size(); ++j)
      cout << jumpList[vacList[i]][j] << " ";
    cout << endl;
  }

  cout << " brutal force: " << endl;
  unordered_map<int, vector<int>> jumpListReference;
  vector<int> vacListReference;
  for (auto&& it = jumpList.begin(); it != jumpList.end(); ++it)
    vacListReference.push_back(it->first);

  for (auto&& i : vacListReference) {
    vector<int> tmpVector;
    for (int j = 0; j < cnf.atoms.size(); ++j) {
      if (cnf.atoms[j].tp == "X")
        continue;
      double dist = calDistPrl(cnf.length, \
                               cnf.atoms[i], \
                               cnf.atoms[j]);
      if (dist <= RCut)
        tmpVector.push_back(j);
    }
    jumpListReference[i] = std::move(tmpVector);
  }
  // for (auto&& i = jumpListReference.begin(); i != jumpListReference.end(); ++i)
  //   vacListReference.push_back(i->first);
  for (int i = 0; i < vacListReference.size(); ++i) {
    cout << vacListReference[i] << " size: " \
         << jumpListReference[vacListReference[i]].size() << endl;
    for (int j = 0; j < jumpListReference[vacListReference[i]].size(); ++j)
      cout << jumpListReference[vacListReference[i]][j] << " ";
    cout << endl;
  }
#endif


}