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

// this class used to extend kmc event to local superbasin event, do not know
// what exact details need to be updated now.
class LSEvent : public KMCEvent {
public:
  LSEvent();
  LSEvent(const pair<int, int>&);
  void exeEvent(Config&, const double&);
};


} // end namespace LS

#endif