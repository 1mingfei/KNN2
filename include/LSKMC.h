#ifndef _LS_KMC_H_
#define _LS_KMC_H_
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
#include "KNUtility.h"
#include "LRUCache.h"

using std::unordered_map;
using std::unordered_set;
using std::vector;
using std::ofstream;
using std::find;

using keras2cpp::Model;
using keras2cpp::Tensor;

class gbCnf;
class LRUCache;

// FPKMC first passage KMC algorithm to speed up
namespace LS {

using arma::vec;
using arma::mat;
typedef vector<double> vd;
typedef vector<vector<double>> vvd;

class LSKMC {
private:
  int& me;
  int& nProcs;
  gbCnf& cnfModifier;
  Config& c0;
  // unordered_map<int, vector<int>>& jumpList;
  unordered_map<string, double>& embedding;

  vector<int>& vacList;
  string& EDiff;
  Model& k2pModelB;
  Model& k2pModelD;
  double& RCut;
  double& RCut2;
  double& temperature;
  double& time;
  double& prefix;
  double& ECutoff;
  double exitTime;
  long long& step;
  ofstream& ofs;
  bool& switchLSKMC;
  bool& switchUnknown;
  vector<string>& elems;
  vector<double>& elemsEffectOffset;
  LRUCache* lru;

  // lists for trapping locations of each atom
  unordered_map<int, unordered_set<int>> trapList;
  // lists for absorbing locations of each atom
  unordered_map<int, unordered_set<int>> absorbList;
  // list of diffusion barriers
  vector<double> barriers;
  // a hashmap trap atom id --> matrix id
  unordered_map<int, int> mapAtomID2MatID;
  unordered_map<int, int> mapMatID2AtomID;

  // event map: i_j --> event_i_j
  unordered_map<string, LSEvent> eventMap;

  vvd VVD_M, VVD_R, VVD_T;
  vd VD_Tau;

  mat Arm_M, Arm_R, Arm_T, Arm_Pi;
  vec Arm_Tau;

  // calculate time of taking each state
  void getVD_Tau(const int&, const int&);

  // helper functions
  // get barrier from hashmap "eventMap" or calculate and put to hashmap
  void getOrPutEvent(const int&, const int&);

  // calculate 2D std vector for M, R, T
  void calVVD_M(const int&);
  void calVVD_R(const int&);
  void calVVD_T(const int&);
  void updateTime();

  // only conduct when trap states are valid
  bool validTrap(const int&);

public:
  LSKMC(int&, \
        int&, \
        gbCnf&, \
        Config&, \
        unordered_map<string, double>&, \
        vector<int>&, \
        string&, \
        Model&, \
        Model&, \
        double&, \
        double&, \
        double&, \
        double&, \
        double&, \
        double& ,\
        long long&, \
        ofstream&, \
        bool&, \
        bool&, \
        vector<string>&, \
        vector<double>&, \
        LRUCache*);

  void testCnfModification();
  static void test_vvd2mat();
  // watch out, this function need to be updated if multiple vacacies is in use
  void searchStatesDFS();

  void helperDFS(const int&, \
                 const int&, \
                 unordered_set<int>&);

  // output surrouding trap states for one vacancy
  void outputTrapCfg(const int&, const string&);
  void outputAbsorbCfg(const int&, const string&);
  void barrierStats();

  void calExitTimePi(const int&);
  int selectEventLSKMC(const int&);
  void executeEvent(const int&, const int&);

};

} // end namespace LS

#endif