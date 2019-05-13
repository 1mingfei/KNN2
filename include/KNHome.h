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
#include <random>
#include <string>
#include <unordered_map>
#include <algorithm> 
#include "gbDef.h"
#include "gbCnf.h"
#include "gbUtl.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::map;
using std::move;
using std::ofstream;
using std::setprecision;
using std::string;
using std::stringstream;
using std::to_string;
using std::unordered_map;
using std::vector;
using std::pair;
using std::min;

class KNHome {
private:
  class gbCnf;

public:
  int me, nProcs;

  unordered_map<string, double> dparams;
  unordered_map<string, int> iparams;
  unordered_map<string, string> sparams;
  unordered_map<string, vector<string>> vsparams;
  unordered_map<string, vector<int>> viparams;



  KNHome(int argc, char* argv[]);
  ~KNHome();

  void parseArgs(int argc, char* argv[]);
  void initParam();
  void readParam();

  void createPreNEB();
};

#endif
