#ifndef _ORDEREDSTRUCT_H_
#define _ORDEREDSTRUCT_H_

#include <vector>
#include <random>

namespace ODS {
using std::pair;
using std::vector;
class OrderedStruct {
public:
  vector<vector<int>> mapping;
  vector<vector<pair<int, int>>> jumpPairs;
public:
  OrderedStruct();
  ~OrderedStruct();
  void generateAuFeOccupInfo();
  void generateCuAuOccupInfo();
  void omit(int, int);
  void makeRandom(int);
  void makeShuffleFraction(int, double);
};
} // end namespace ODS

#endif //_ORDEREDSTRUCT_H_
