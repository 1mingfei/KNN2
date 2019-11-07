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
#include <set>
#include <algorithm> 

// #include "KNHome.h"
#include "gbDef.h"

using std::pair;


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
  void calRate(const Config&, const double&);
  void calProb(const double&);
  void setcProb(const double&);
  void exeEvent(Config&);
// friend class KNHome;

};
#endif