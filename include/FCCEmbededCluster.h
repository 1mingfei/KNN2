#ifndef _FCC_EMBEDED_CLUSTER_H_
#define _FCC_EMBEDED_CLUSTER_H_

#include "KNHome.h"

namespace FCCEmbededCluster {

struct occupInfo_256 {
  vector<vector<int>> mapping;
  vector<vector<pair<int, int>>> jumpPairs;

  occupInfo_256() {
    mapping = vector<vector<int>> (7, vector<int> (256, 0) );
    // possibleJumpPairs = vector<vector<pair<int, int>>> (7, \
    //                            vector<pair<int, int>>);

    // L10
    for (int i : {0,1,4,5,16,17,20,21,64,68,65,69,80,84,81,85})
      mapping[0][i] = 1;
    for (int i : {3,19,2,18,7,23,6,22,67,83,66,82,71,87,70,86})
      mapping[0][i] = 2;
    jumpPairs.emplace_back(vector<pair<int, int>> ({ {82, 83},
                                                     {82, 80},
                                                     {82, 148},
                                                     {83, 160},
                                                     {83, 81},
                                                     {18, 23},
                                                     {18, 84} }) );


  };
  ~occupInfo_256() {};
};

} // namespace FCCEmbededCluster

#endif