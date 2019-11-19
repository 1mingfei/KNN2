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
                        unordered_map<int, \
                        vector<int>>& jumpList,\
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
  /* 3 */
  std::swap(cnf.atoms[first].NBL, cnf.atoms[second].NBL);
  /* 4 */
  for (int i : cnf.atoms[first].NBL) {
    if (find(cnf.atoms[first].NBL.begin(), \
             cnf.atoms[first].NBL.end(), \
             second) == cnf.atoms[first].NBL.end()) {
      replace(cnf.atoms[i].NBL.begin(), cnf.atoms[i].NBL.end(), \
              first, second);
    }
  }
  for (int i : cnf.atoms[second].NBL) {
    if (find(cnf.atoms[second].NBL.begin(), \
             cnf.atoms[second].NBL.end(), \
             first) == cnf.atoms[second].NBL.end()) {
    replace(cnf.atoms[i].NBL.begin(), cnf.atoms[i].NBL.end(), \
            second, first);
    } 
  }
  /* 5: to do update jumpList by recalculating atoms within cutoff of 3.0A */
  // find the previous vac; update its neighbor list
  vector<int> res;
  for (int i = 0; i < cnf.atoms.size(); ++i) {
    double dist = calDistPrl(cnf.length, \
                             cnf.atoms[first], \
                             cnf.atoms[i]);
    // if ((dist <= RCut) && (i % cnf.atoms.size() != first))
    if ((dist <= RCut) && (i != first))
      res.push_back(i);
  }

  jumpList[first] = std::move(res);

}