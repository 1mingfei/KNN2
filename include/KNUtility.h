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

typedef vector<vector<double>> vvd;

mat vvd2mat(vector<vector<double>>& vIn);

#endif