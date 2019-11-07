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
  pair<int, int> jumpPair;

public:
  KMCEvent();
  KMCEvent(const pair<int, int>&);
  ~KMCEvent();
  void getRate(const Config&);
  void calRate(const double&, const double&);
  void exeEvent(Config&);
// friend class KNHome;

};
#endif