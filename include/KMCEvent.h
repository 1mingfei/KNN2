#ifndef _KMC_EVENT_H_
#define _KMC_EVENT_H_

#include "KNHome.h"

class KNHome::KMCEvent {
  KNHome& kn;
  unordered_map<string, double>& dparams;
  unordered_map<string, int>& iparams;
  unordered_map<string, string>& sparams;
  unordered_map<string, vector<string>>& vsparams;
  unordered_map<string, vector<int>>& viparams;

public:
  KMCEvent(KNHome&);
  ~KMCEvent();



};
#endif