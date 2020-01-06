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

class KNHome;

class KMCEvent {
private:
  double barrier;
  double rate;
  double energyChange;
  double prob;
  double cProb; // cumulative probability
  pair<int, int> jumpPair;

public:
  KMCEvent();
  KMCEvent(const pair<int, int>&);
  ~KMCEvent();

  double getBarrier() const;
  double getRate() const;
  double getProb() const;
  double getcProb() const;
  double getEnergyChange() const;
  pair<int, int> getJumpPair() const;


  // void calRate(const Config&, const double&, const double&);
  void setBarrier(const double&);
  void setRate(const double&);
  void calProb(const double&);
  void setcProb(const double&);
  void setEnergyChange(const double&);
  void exeEvent(Config&, const double&);

};
#endif