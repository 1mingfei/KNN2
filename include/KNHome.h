#ifndef _KNHome_H_
#define _KNHome_H_

#include <math.h>
#include <mpi.h>
// #include <omp.h>

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
#include <list>
#include <algorithm>
#include "armadillo"
#include "gbDef.h"
#include "KMCEvent.h"
#include "LSEvent.h"
#include "gbUtl.h"
#include "model.h"
#include "LSKMC.h"
#include "KNUtility.h"
#include "OrderedStruct.h"
#include "LRUCache.h"

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
using std::list;
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
class LRUCache;

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
  bool switchLSKMC, switchUnknown;
  unordered_map<string, double> embedding;
  unordered_map<string, int> eventListMap;
  vector<string> elems;
  vector<double> elemsEffectOffset;

  ofstream ofs;

  string EDiff;

  bool isTrapped(const double&);
  double updateTime();

public:
  int me, nProcs;

  Model k2pModelB;
  Model k2pModelD;

  unordered_map<string, double> dparams;
  unordered_map<string, long long> iparams;
  unordered_map<string, bool> bparams;
  unordered_map<string, string> sparams;
  unordered_map<string, vector<string>> vsparams;
  unordered_map<string, vector<int>> viparams;
  unordered_map<string, vector<double>> vdparams;

  long long LRUSize;
  LRUCache* lru;

  KNHome(int argc, char* argv[]);
  ~KNHome();

  /* KNParam.cpp */
  void parseArgs(int argc, char* argv[]);
  void initParam();
  void readParam();

  /* KNSolidSol.cpp */
  void createPreNEB(gbCnf&);

  /* KNRandom.cpp */
  void createRandom(gbCnf&, \
                    const int&, \
                    const int&,  \
                    const string&, \
                    const double&, \
                    const vector<int>&, \
                    const vector<int>&);

  void createRandomUniform(gbCnf&, \
                           const int&, \
                           const int&,  \
                           const string&, \
                           const double&, \
                           const vector<int>&, \
                           const vector<int>&);

  void createRandomSpecific(gbCnf&, \
                            const int&, \
                            const int&,  \
                            const string&, \
                            const double&, \
                            const vector<int>&, \
                            const vector<int>&);

  /* KNOrdered.cpp */
  int createSingle(const int&, \
                   int, \
                   gbCnf&, \
                   const vector<int>&, \
                   const double&, \
                   const string&, \
                   const ODS::OrderedStruct&, \
                   const pair<string, string>&);

  void createOrdered(gbCnf&, \
                     const vector<int>&, \
                     const double&, \
                     const string&);

  void createOrderedRandom(gbCnf&, \
                           const vector<int>&, \
                           const double&, \
                           const string&, \
                           const int&);

  void createOrderedDiffCon(gbCnf&, \
                            const vector<int>&, \
                            const double&, \
                            const string&, \
                            const int&);

  void createOrderedAntiPhase(gbCnf&, \
                              const vector<int>&, \
                              const double&, \
                              const string&, \
                              const int&);
  /* KNvasp.cpp */
  void prepVASPFiles(const string&, \
                     const vector<int>&, \
                     const map<string, int>&, \
                     const string&);

  /* KNEncode.cpp */
  void KNEncode(gbCnf&);
  vector<vector<int>> readPairs(const string&);

  /* KNBondCount.cpp */
  void KNBondCount(gbCnf&);

  /* KMCSimulation.cpp */
  void KMCInit(gbCnf&);
  void getVacList();
  void KMCSimulation(gbCnf&);
  void buildEmbedding();
  vector<double> calRate(Config&, \
                         const double&, \
                         gbCnf&, \
                         pair<int, int>);

  void buildEventList_serial(gbCnf&);
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
  void findClts(gbCnf&, const string&);
  void loopConfig(gbCnf&);

};


// random generator function:
inline int myRandInt(int minVal, int maxVal) {
  return minVal + (std::rand() % static_cast<int>(maxVal - minVal + 1));
}

#endif