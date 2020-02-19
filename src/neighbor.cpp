#include "gbCnf.h"
#include "KNUtility.h"
#define FNN_DIST 3.0

void gbCnf::getNBL_serial(Config& cnf, double Rcut = 3.8) {
  vector<double> tmpLength;
  tmpLength = cnf.length;
  vector<KNAtom> tmpAtoms = cnf.atoms;
  for (int i = 0; i < cnf.atoms.size(); ++i) {
    vector<int> res;
    int k = 0;
    for (int j = 0; j < tmpAtoms.size(); ++j) {
      double dist = calDistPrl(tmpLength, tmpAtoms[i], tmpAtoms[j]);
      if ((dist <= Rcut) && (j % cnf.atoms.size() - i != 0)) {
        res.push_back(j);
        if (dist <= FNN_DIST)
          cnf.atoms[i].FNNL[k++] = j;
      }
    }
    cnf.atoms[i].NBL = std::move(res);
  }
}

void gbCnf::getNBL(Config& cnf, double Rcut = 3.8) {
  vector<double> tmpLength;
  tmpLength = cnf.length;
  vector<KNAtom> tmpAtoms = cnf.atoms;

  int size = cnf.atoms.size();
  int quotient = size / nProcs;
  int remainder = size % nProcs;
  int nCycle = remainder ? (quotient + 1) : quotient;


  // allocate largest memory for all nAtom (this is slightly larger than needed)
  // Note that "nCycle * nProcs >= nAtoms"
  // all data store here after this
  int** NBLArry;
  int** FNNLArry;

  NBLArry = alloc_2d_array<int>(nCycle * nProcs, 18);
  FNNLArry = alloc_2d_array<int>(nCycle * nProcs, 12);

  // smallest buff for gathering in each cycle
  int* buff_NBLArry = new int [18];
  int* buff_FNNLArry = new int [12];


  for (int j = 0; j < nCycle; ++j) {

    for (int i = (j * nProcs); i < ((j + 1) * nProcs); ++i) {

      if ((me == 0) && (i % nProcs != 0)) {
        MPI_Recv(&NBLArry[i][0], 18, MPI_INT, (i % nProcs), 0, \
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&FNNLArry[i][0], 12, MPI_INT, (i % nProcs), 1, \
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }

      if (i % nProcs != me) continue;

      for (int l = 0; l < 18; ++l) {
        buff_NBLArry[l] = -1;
        if (l >= 12)
          continue;
        buff_FNNLArry[l] = -1;
      }

      if (i >= size) continue;

      int num_NBL = 0, num_FNNL = 0;
      for (int l = 0; l < tmpAtoms.size(); ++l) {
        double dist = calDistPrl(tmpLength, tmpAtoms[i], tmpAtoms[l]);
        if ((dist <= Rcut) && (l % size - i != 0)) {
          buff_NBLArry[num_NBL++] = l;
          if (dist <= FNN_DIST)
            buff_FNNLArry[num_FNNL++] = l;
        }
      }

      if (me != 0) {
        MPI_Send(&buff_NBLArry[0], 18, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&buff_FNNLArry[0], 12, MPI_INT, 0, 1, MPI_COMM_WORLD);
      } else {
        for (int ii = 0; ii < 18; ++ii) {
          NBLArry[j * nProcs + i % nProcs][ii] = buff_NBLArry[ii];
          if (ii >= 12) continue;
          FNNLArry[j * nProcs + i % nProcs][ii] = buff_FNNLArry[ii];
        }
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  MPI_Bcast(&NBLArry[0][0], nCycle * nProcs * 18, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&FNNLArry[0][0], nCycle * nProcs * 12, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < 18; ++j) {
      cnf.atoms[i].NBL.push_back(NBLArry[i][j]);
      if (j >= 12) continue;
      cnf.atoms[i].FNNL[j] = FNNLArry[i][j];
    }
  }

  // for (int j = 0; j < nProcs; ++j) {
  //   if (me == j) {
  //     cout << "proc #" << j << " " << cnf.atoms.size() << endl;
  //     for (int i = 499 ; i > 490; --i) {
  //       cout << "proc #" << j << " atom #" << i << " NB size: "
  //            << cnf.atoms[i].NBL.size() << "\n";
  //     }

  //     for (int i = 499 ; i > 490; --i) {
  //       cout << "proc #" << j << "\n";
  //       for (int k = 0; k < 12; ++k)
  //         cout << " " << cnf.atoms[i].FNNL[k] << " ";
  //       cout << "\n";
  //     }
  //   }
  // }

  delete [] buff_NBLArry;
  delete [] buff_FNNLArry;

  free(NBLArry[0]);
  free(FNNLArry[0]);
  free(NBLArry);
  free(FNNLArry);
}

int gbCnf::getExpdParam(const Config& cnf, const double Rcut = 3.8) {
  if (cnf.length[Z] > 2.0 * Rcut) {
    return 1;
  } else {
    return(static_cast<int>((2.0 * Rcut / cnf.length[Z]) + 1));
  }
}

/* expand cell in +/- Z direction */
vector<KNAtom> gbCnf::expandCellZ(const Config& cnf, const int factor) {
  vector<KNAtom> res;
  int initSize = cnf.atoms.size();
  for (int i = 0; i < initSize; ++i) {
    res.push_back(cnf.atoms[i]);
  }
  for (int i = initSize; i < initSize * factor; ++i) {
    KNAtom atm = cnf.atoms[i % initSize];
    res.push_back(KNAtom(i, atm.tp, atm.pst[X], atm.pst[Y],
           atm.pst[Z] + cnf.length[Z]*int(i / initSize)));
  }
  return res;
}

/*calculate distance between one atom in configuration and one from ref*/
double gbCnf::calDist(const vector<double> length, \
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
double gbCnf::calDistPrl(const vector<double>& length, \
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

  double dist = sqrt(a * a + b * b + c * c);
  return dist;
}
