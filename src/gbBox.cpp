#include "gbCnf.h"
#include "KNHome.h"

#define DOMLEN 4.34
#define D_SIGMA 0.15
#define D_CUTOFF 0.4

void gbCnf::initBox(Config& c) {
  crossProd33(c.bvy, c.bvz, c.tvx);
  crossProd33(c.bvz, c.bvx, c.tvy);
  crossProd33(c.bvx, c.bvy, c.tvz);

  double vol = vecInnProd33(c.bvx, c.tvx);
  double inv = 1. / vol;

  // normalize
  scaleVec(c.tvx, inv);
  scaleVec(c.tvy, inv);
  scaleVec(c.tvz, inv);

  double iheight[3];
  iheight[0] = std::sqrt(square33(c.tvx));
  iheight[1] = std::sqrt(square33(c.tvy));
  iheight[2] = std::sqrt(square33(c.tvz));
  // for (auto&& i : {X, Y, Z})
  //   c.scale[i] = static_cast<int>(std::ceil(rcut * iheight[i]));
}

void gbCnf::cnvprl2pst(Config& c) {
  for (KNAtom& atm : c.atoms) {
    for (const int i : {0, 1, 2})
      atm.pst[i] =
          atm.prl[0] * c.bvx[i] + atm.prl[1] * c.bvy[i] + atm.prl[2] * c.bvz[i];
  }
}

void gbCnf::wrapAtomPos(Config& tmpc) {

  vector<double> boxlo(3), boxhi(3), npst(3);
  for (const int& k : {0, 1, 2}) {
    boxlo[k] = 0.0 * tmpc.bvx[k] + 0.0 * tmpc.bvy[k] + 0.0 * tmpc.bvz[k];
    boxhi[k] = 1.0 * tmpc.bvx[k] + 1.0 * tmpc.bvy[k] + 1.0 * tmpc.bvz[k];
  }

  for (auto&& atm : tmpc.atoms) {
    for (int ix = -1; ix <= 1; ix++) {
      for (int iy = -1; iy <= 1; iy++) {
        for (int iz = -1; iz <= 1; iz++) {
          if ((ix == 0) && (iy == 0) && (iz == 0)) continue;
          npst[0] = atm.pst[0] + ix * tmpc.bvx[0] + iy * tmpc.bvy[0] +
                    iz * tmpc.bvz[0];
          if (npst[0] > boxhi[0] || npst[0] < boxlo[0]) continue;

          npst[1] = atm.pst[1] + ix * tmpc.bvx[1] + iy * tmpc.bvy[1] +
                    iz * tmpc.bvz[1];
          if (npst[1] > boxhi[1] || npst[1] < boxlo[1]) continue;

          npst[2] = atm.pst[2] + ix * tmpc.bvx[2] + iy * tmpc.bvy[2] +
                    iz * tmpc.bvz[2];
          if (npst[2] > boxhi[2] || npst[2] < boxlo[2]) continue;
          // update
          for (const int& it : {0, 1, 2}) atm.pst[it] = npst[it];
        }
      }
    }
  }
}
/*convert cart pos to frac pos*/
void gbCnf::cnvpst2prl(Config& c) {
  /*
  for (KNAtom& atm : c.atoms) {
    for (const int i : {0, 1, 2})
      atm.pst[i] =
          atm.prl[0] * c.bvx[i] + atm.prl[1] * c.bvy[i] + atm.prl[2] * c.bvz[i];
  }
  */
  mat A = {{c.bvx[0], c.bvx[1], c.bvx[2]},
           {c.bvy[0], c.bvy[1], c.bvy[2]},
           {c.bvz[0], c.bvz[1], c.bvz[2]}};

  for (KNAtom& atm : c.atoms) {
    vec X;
    vec B = {atm.pst[0],
             atm.pst[1],
             atm.pst[2]};

    //X = solve(A, B, arma::solve_opts::allow_ugly);
    X = solve(A, B);

    atm.prl[0] = X[0];
    atm.prl[1] = X[1];
    atm.prl[2] = X[2];
  }
}

void gbCnf::wrapAtomPrl(Config& tmpc) {
  for (auto&& atm : tmpc.atoms) {
    for (int i = 0; i < 3; ++i) {
      double factor = std::floor(atm.prl[i]);
      atm.prl[i] -= factor;
    }
  }
}

void gbCnf::perturb(Config& cnf) {

  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0.0, D_SIGMA);

  for (auto&& atm : cnf.atoms) {
    for (const int i : {0, 1, 2}) {
      double displacement = distribution(generator);
      while (displacement > D_CUTOFF) {
        displacement = distribution(generator);
      }
      atm.pst[i] += displacement;
    }
  }
  cnvpst2prl(cnf);
}