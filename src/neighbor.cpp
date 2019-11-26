/*
 * Author: 1mingfei 
 * Date:   2019-04-17
 * Purpose: 
 */

#include "gbCnf.h"

void KNHome::gbCnf::getNBL(Config& cnf, double Rcut = 3.8) {
  //int factor = getExpdParam(cnf, Rcut);
  vector<double> tmpLength;
  tmpLength = cnf.length;
  //tmpLength[Z] *= double(factor);
  //vector<KNAtom> tmpAtoms = expandCellZ(cnf, factor);
  
  vector<KNAtom> tmpAtoms = cnf.atoms;
  for (int i = 0; i < cnf.atoms.size(); ++i) {
    vector<int> res;
    for (int j = 0; j < tmpAtoms.size(); ++j) {
      double dist = calDistPrl(tmpLength, tmpAtoms[i], tmpAtoms[j]);
      if ((dist <= Rcut) && (j % cnf.atoms.size() - i != 0))
        res.push_back(j);
    }
    cnf.atoms[i].NBL = std::move(res);
  }
}

int KNHome::gbCnf::getExpdParam(const Config& cnf, const double Rcut = 3.8) {
  if (cnf.length[Z] > 2.0*Rcut) {
    return 1;
  } else {
    return(static_cast<int>((2.0*Rcut/cnf.length[Z])+1));
  }
}

/* expand cell in +/- Z direction */
vector<KNAtom> KNHome::gbCnf::expandCellZ(const Config& cnf, const int factor) {
  vector<KNAtom> res;
  int initSize = cnf.atoms.size();
  for (int i = 0; i < initSize; ++i) {
    res.push_back(cnf.atoms[i]);
  }
  for (int i = initSize; i < initSize*factor; ++i) {
    KNAtom atm = cnf.atoms[i%initSize];
    res.push_back(KNAtom(i, atm.tp, atm.pst[X], atm.pst[Y],
           atm.pst[Z] + cnf.length[Z]*int(i/initSize)));
  }
  return res;
}

/*calculate distance between one atom in configuration and one from ref*/
double KNHome::gbCnf::calDist(const vector<double> length, \
                              const KNAtom& atm1, \
                              const KNAtom& atm2) {
  double xi = atm1.pst[X];
  double xj = atm2.pst[X];
  double yi = atm1.pst[Y];
  double yj = atm2.pst[Y];
  double zi = atm1.pst[Z];
  double zj = atm2.pst[Z];
  double a, b, c;

  if (xj - xi >= 0.5 * length[X]) {
    a = (xi - xj + length[X]);
  } else if (xj - xi <  -0.5 * length[X]) {
    a = (xi - xj - length[X]); 
  } else {
    a = xi - xj;
  }

  b = yi - yj;

  if (zj - zi >= 0.5 * length[Z]) {
    c = (zi - zj + length[Z]);
  } else if (zj - zi <  -0.5 * length[Z]) {
    c = (zi - zj - length[Z]); 
  } else {
    c = zi - zj;
  }

  double dist = sqrt(a*a + b*b + c*c);
  return dist;
}

/*calculate distance between one atom in configuration and one from ref*/
double KNHome::gbCnf::calDistPrl(const vector<double>& length, \
                                 const KNAtom& atm1, \
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
