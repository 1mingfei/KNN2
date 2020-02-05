#ifndef _KNHome_H_
#define _KNHome_H_

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
#include <set>
#include <array>
#include <vector>
#include <queue>
#include <algorithm>
#include "armadillo"
#include "gbDef.h"
// #include "gbCnf.h"
#include "KMCEvent.h"
#include "LSEvent.h"
#include "gbUtl.h"
#include "model.h"
#include "LSKMC.h"
#include "KNUtility.h"
#include "FCCEmbededCluster.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::map;
using std::multimap;
using std::unordered_multimap;
using std::queue;
using std::move;
using std::ofstream;
using std::setprecision;
using std::string;
using std::stringstream;
using std::to_string;
using std::unordered_map;
using std::vector;
using std::set;
using std::pair;
using std::min;
using std::distance;
using std::replace;
using std::swap;
using std::make_pair;

using arma::vec;
using arma::mat;

using keras2cpp::Model;
using keras2cpp::Tensor;

class gbCnf;

class KNHome {
private:

  Config c0;
  vector<KMCEvent> eventList;
  vector<int> vacList;
  // unordered_map<int, vector<int>> jumpList;
  // unordered_map<int, vector<int>> oldJumpList;

  double RCut, RCut2;
  double temperature;
  double time;
  double prefix;
  double kTot;
  double E_tot; // total energy change of the system
  double ECutoff;
  long long maxIter, iter;
  long long step, trapStep;
  int nTallyConf, nTallyOutput;
  unordered_map<string, double> embedding;
  unordered_map<string, int> eventListMap;

  bool switchEngy;

  bool isTrapped(const double&);
  double updateTime();

public:
  int me, nProcs;

  Model k2pModelB;
  Model k2pModelD;

  unordered_map<string, double> dparams;
  unordered_map<string, int> iparams;
  unordered_map<string, bool> bparams;
  unordered_map<string, string> sparams;
  unordered_map<string, vector<string>> vsparams;
  unordered_map<string, vector<int>> viparams;


  KNHome(int argc, char* argv[]);
  ~KNHome();

  /* KNParam.cpp */
  void parseArgs(int argc, char* argv[]);
  void initParam();
  void readParam();

  /* KNSolidSol.cpp */
  void createPreNEB();

  /* KNOrdered.cpp */
  void createOrdered(gbCnf&, \
                     const vector<int>&, \
                     const double&, \
                     const string&);

  int createOrderedSingle(const int&, \
                           int, \
                           gbCnf&, \
                           const vector<int>&, \
                           const double&, \
                           const string&, \
                           const FCCEmbededCluster::occupInfo_256&);

  /* KNvasp.cpp */
  void prepVASPFiles(const string&, \
                     const vector<int>&, \
                     const map<string, int>&, \
                     const string&);

  /* KNEncode.cpp */
  void KNEncode();
  vector<vector<int>> readPairs(const string&);

  /* KNBondCount.cpp */
  void KNBondCount();

  /* KMCSimulation.cpp */
  void KMCInit(gbCnf&);
  void getVacList();
  void KMCSimulation(gbCnf&);
  void buildEmbedding();
  vector<double> calRate(Config&, \
                         const double&, \
                         gbCnf&, \
                         pair<int, int>);

  void buildEventList(gbCnf&);
  void updateEventList(gbCnf&, \
                       const pair<int, int>&, \
                       const int&);

  KMCEvent selectEvent(int&);
  void updateEnergy(const int&);

  /* test keras2cpp */
  void testK2P();

  /* LSKMC */
  void LSKMCOneRun(gbCnf&);
  void LSKMCSimulation(gbCnf&);

  /* findClusters */
  void findClts(gbCnf&);
};

#endif
