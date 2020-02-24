//
// Created by Zhucong Xi on 2/19/20.
//

#include "OrderedStruct.h"

OrderedStruct::OrderedStruct() {
}
OrderedStruct::~OrderedStruct() {
}
void OrderedStruct::generateAuFeOccupInfo() {
  mapping = vector<vector<int>>(6, vector<int>(256, 0));

  // L10
  for (int i : {0, 1, 4, 5, 16, 17, 20, 21, 64, 68, 65, 69, 80, 84, 81, 85})
    mapping[0][i] = 1;
  for (int i : {3, 19, 2, 18, 7, 23, 6, 22, 67, 83, 66, 82, 71, 87, 70, 86})
    mapping[0][i] = 2;
  jumpPairs.emplace_back(vector<pair<int, int>>({{82, 83},
                                                 {82, 80},
                                                 {82, 148},
                                                 {83, 160},
                                                 {83, 81},
                                                 {83, 144},
                                                 {18, 23},
                                                 {18, 84},
                                                 {18, 19}
                                                }));


  // L12
  for (int i : {0, 4, 16, 20, 64, 68, 80, 84})
    mapping[1][i] = 1;
  for (int i : {1, 3, 5, 17, 19, 2, 18, 7, 23, 6, 21, 22, 67, 83, 65, 66, 69,
                81, 82, 85, 71, 87, 70, 86})
    mapping[1][i] = 2;
  jumpPairs.emplace_back(vector<pair<int, int>>({{82, 83},
                                                 {82, 80},
                                                 {82, 148},
                                                 {83, 160},
                                                 {83, 81},
                                                 {83, 144},
                                                 {18, 23},
                                                 {18, 84},
                                                 {18, 19}
                                                }));

  // L10*
  for (int i : {0, 2, 4, 5, 16, 18, 20, 21, 64, 68, 66, 69, 80, 84, 82, 85})
    mapping[2][i] = 1;
  for (int i : {3, 19, 1, 17, 7, 23, 6, 22, 67, 83, 65, 81, 71, 87, 70, 86})
    mapping[2][i] = 2;
  jumpPairs.emplace_back(vector<pair<int, int>>({{82, 83},
                                                 {82, 80},
                                                 {82, 148},
                                                 {83, 160},
                                                 {83, 81},
                                                 {83, 144},
                                                 {18, 23},
                                                 {18, 84},
                                                 {18, 19},
                                                 {1, 4},
                                                 {1, 2},
                                                 {1, 194},
                                                 {66, 65},
                                                 {66, 71},
                                                 {66, 129}
                                                }));

  // W2
  // purple
  for (int i : {8, 16, 17, 32, 36,
      // layer 2
                3, 18, 23, 38, 2, 7, 22, 27,
      // layer 3
                68, 69, 88, 72, 96,
      // layer 4
                83, 98, 67, 82, 87, 102,
      // layer 5
                128, 129, 148, 149, 168, 132, 133, 152
  })
    mapping[3][i] = 1;
  // yellow
  for (int i : {0, 1, 20, 21, 40, 4, 5, 24,
      // layer 2
                11, 6, 19, 34,
      // layer 3
                64, 65, 84, 85, 104, 100, 81, 80,
      //layer 4
                66, 71, 86, 91, 70, 75,
      //layer 5
                136, 144, 145, 164, 160
  })
    mapping[3][i] = 2;
  jumpPairs.emplace_back(vector<pair<int, int>>({
                                                    {84, 65},
                                                    {84, 87},
                                                    {148, 151},
                                                    {148, 135},
                                                    {148, 87},
                                                    {148, 71},
                                                    {148, 129},
                                                    {148, 145},
                                                    {148, 133},
                                                    {20, 18},
                                                    {20, 215},
                                                    {20, 21},
                                                    {20, 17},
                                                    {20, 5},
                                                    {75, 72},
                                                    {75, 152},
                                                    {75, 74},
                                                    {75, 69},
                                                    {75, 70}
                                                }));

  // Z1
  for (int i : {0, 1, 4, 5, 16, 17, 20, 21})
    mapping[4][i] = 1;
  for (int i : {3, 19, 2, 18, 7, 23, 6, 22, 67, 83, 66, 82, 71, 87, 70, 86,
                64, 68, 65, 69, 80, 84, 81, 85})
    mapping[4][i] = 2;
  jumpPairs.emplace_back(vector<pair<int, int>>({{82, 83},
                                                 {82, 80},
                                                 {82, 148},
                                                 {83, 160},
                                                 {83, 81},
                                                 {83, 144},
                                                 {18, 23},
                                                 {18, 84},
                                                 {18, 19}
                                                }));
  // L12*
  for (int i : {0, 16, 64, 80, 7, 11, 23, 27, 12, 28, 76, 92})
    mapping[5][i] = 1;
  for (int i : {3, 19, 67, 83, 17, 1, 2, 18, 65, 81, 66, 82, 4, 20, 68, 84,
                87, 71, 5, 21, 6,
                22, 69, 70, 85, 86, 8, 24, 72, 88, 75, 91, 9, 10, 25, 26, 73,
                74, 89, 90, 15,
                31, 79, 95, 94, 93, 78, 77, 30, 29, 14, 13})
    mapping[5][i] = 2;
  jumpPairs.emplace_back(vector<pair<int, int>>({{82, 83},
                                                 {82, 80},
                                                 {82, 65},
                                                 {82, 84},
                                                 {82, 148},
                                                 {83, 160},
                                                 {83, 81},
                                                 {83, 80},
                                                 {79, 137},
                                                 {79, 76},
                                                 {79, 141},
                                                 {79, 74}
                                                }));
}

void OrderedStruct::generateCuAuOccupInfo() {
  mapping = vector<vector<int>>(1, vector<int>(256, 0));
  // Mg 1 or Al 0 as open circles;
  // Zn 2 atoms as filled circles,
  // disorder at these sites 0, 1, 2 as hatched circles

  // open circles
  for (int i : {2, 3, 18, 19, 66, 67, 82, 83, 5, 8, 9, 21, 24, 25, 69, 72, 73,
                85, 88, 89, 14, 30, 78, 94})
    mapping[0][i] = rand() % 2;

  // filled circles
  for (int i : {0, 1, 64, 65, 16, 17, 80, 81, 6, 11, 10, 22, 27, 26, 70, 75,
                74, 86, 91, 90, 13, 29, 77, 93})
    mapping[0][i] = 2;

  // hatched circles
  for (int i : {4, 20, 7, 23, 68, 84, 71, 87, 79, 76, 15, 12, 95, 92, 31, 28})
    mapping[0][i] = rand() % 3;

  jumpPairs.emplace_back(vector<pair<int, int>>({{71, 132},
                                                 {71, 66},
                                                 {71, 70},
                                                 {71, 84},
                                                 {71, 65},
                                                 {22, 11},
                                                 {22, 69},
                                                 {22, 23},
                                                 {18, 23},
                                                 {18, 84},
                                                 {18, 19},
                                                 {20, 18},
                                                 {20, 215},
                                                 {20, 21},
                                                 {20, 17},
                                                 {20, 5}
                                                }));
}
// change purple or yellow to base element, 1 for purple, 2 for yellow
// i is the sturcture index
void OrderedStruct::omit(int i, int colorI) {
  for (int j = 0; j < mapping[i].size(); ++j) {
    if (mapping[i][j] == colorI) {
      mapping[i][j] = 0;
    }
  }
}
// This function re-partitions yellow and purple atoms randomly.
// i is the structure index
void OrderedStruct::makeRandom(int i) {
  for (int& pos : mapping[i]) {
    if (pos != 0) {
      pos = (rand() % 2) + 1;
    }
  }
}
void OrderedStruct::makeShuffleFraction(int i,
                                        double purpleFraction) {
  vector<int> indexMap;
  for (int j = 0; j < mapping[i].size(); ++j) {
    if (mapping[i][j] != 0) {
      mapping[i][j] = 2;
      indexMap.push_back(j);
    }
  }
  auto purpleNum = static_cast<int>(round(purpleFraction * indexMap.size()));
  shuffle(indexMap.begin(),
          indexMap.end(),
          std::default_random_engine(rand()));
  for (int k = 0; k < purpleNum; ++k) {
    mapping[i][indexMap[k]] = 1;
  }
}
