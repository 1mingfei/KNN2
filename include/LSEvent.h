#ifndef _LS_EVENT_H_
#define _LS_EVENT_H_

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
#include "KMCEvent.h"


using std::pair;
using std::unordered_map;
using std::unordered_set;

namespace LS {

class LSEvent : public KMCEvent {
public:
  LSEvent(const pair<int, int>&);
  void exeEvent(Config&, const double&);
};


} // end namespace LS

#endif