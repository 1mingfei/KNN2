/*
 * Author: 1mingfei 
 * Date:   2019-05-09
 * Purpose: function to generate FCC structures 
 * self-explained
 */

#include "gbCnf.h"
Config KNHome::gbCnf::getFCCConv(const double LC, const string elem, 
                                 const vector<int>& factors) {
  Config cnf;
  cnf.length[X] = cnf.bvx[X] = LC * factors[X];
  cnf.length[Y] = cnf.bvy[Y] = LC * factors[Y];
  cnf.length[Z] = cnf.bvz[Z] = LC * factors[Z];
  int c = 0;
  double na = static_cast<double>(factors[X]);
  double nb = static_cast<double>(factors[Y]);
  double nc = static_cast<double>(factors[Z]);
  for (int k = 0; k < factors[Z]; ++k) {
    for (int j = 0; j < factors[Y]; ++j) {
      for (int i = 0; i < factors[X]; ++i) {
        double xref = static_cast<double>(i);
        double yref = static_cast<double>(j);
        double zref = static_cast<double>(k);
        double x, y, z;
        for (int l = 0; l < 4; ++l) {
          switch (l % 4) {
            case 1 : x = xref + 0.5;
                     y = yref + 0.5;
                     z = zref;
                     break;
            case 2 : x = xref + 0.5;
                     y = yref;
                     z = zref + 0.5;
                     break;
            case 3 : x = xref;
                     y = yref + 0.5;
                     z = zref + 0.5;
                     break;
            default : x = xref;
                      y = yref;
                      z = zref;
                      break;
          }
          x /= na;
          y /= nb;
          z /= nc;
          KNAtom atm(c, elem, x, y, z);
          cnf.atoms.push_back(atm);
          ++c;
        }
      }
    }
  }
  cnf.natoms = c;
  return cnf;
}
