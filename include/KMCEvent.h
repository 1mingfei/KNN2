#ifndef _KMC_EVENT_H_
#define _KMC_EVENT_H_

#include <math.h>
#include <mpi.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <algorithm>

// #include "KNHome.h"
#include "gbDef.h"


using std::pair;
using std::unordered_map;
using std::unordered_set;
using std::find;
using std::remove;

class KNHome;

class KMCEvent {
private:
  double rate;
  double prob;
  double cProb; // cumulative probability
  pair<int, int> jumpPair;

public:
  KMCEvent();
  KMCEvent(const pair<int, int>&);
  ~KMCEvent();
  double getRate() const;
  double getProb() const;
  double getcProb() const;
  pair<int, int> getJumpPair() const;

  // void calRate(const Config&, const double&, const double&);
  void setRate(const double&);
  void calProb(const double&);
  void setcProb(const double&);
  void exeEvent(Config&, unordered_map<int, vector<int>>&, const double&);
// friend class KNHome;

};
#endif