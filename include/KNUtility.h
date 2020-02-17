#ifndef _KN_UTILITY_H_
#define _KN_UTILITY_H_

#include <math.h>
#include <mpi.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <algorithm>
// #include <bits/stdc++.h>
#include "armadillo"
#include "gbDef.h"
#include "gbCnf.h"
#include "KNHome.h"

using std::unordered_map;
using std::unordered_set;
using std::vector;
using std::ofstream;
using arma::mat;
using arma::vec;

typedef vector<double> vd;
typedef vector<vector<double>> vvd;

mat vvd2mat(vvd&);
vec vd2vec(vd&);
vd mat2vd(mat&);

template <class T>
T** alloc_2d_array(int rows, int cols) {
  T* data = (T*) malloc(rows * cols * sizeof(T));
  T** array= (T**) malloc(rows * sizeof(T*));
  for (int i = 0; i < rows; ++i)
      array[i] = &(data[cols * i]);
  return array;
}

#endif