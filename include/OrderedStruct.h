//
// Created by Zhucong Xi on 2/19/20.
//

#ifndef _ORDEREDSTRUCT_H_
#define _ORDEREDSTRUCT_H_

#include "KNHome.h"

namespace FCCCluster {
class OrderedStruct {
 public:
  vector<vector<int>> mapping;
  vector<vector<pair<int, int>>> jumpPairs;
 public:
  OrderedStruct();
  virtual ~OrderedStruct();
  void generateAuFeOccupInfo();
  void generateCuAuOccupInfo();
  void omit(int, int);
  void makeRandom(int);
  void makeShuffleFraction(int, double);
};
}
#endif //_ORDEREDSTRUCT_H_
