#include "KNUtility.h"

mat vvd2mat(vector<vector<double>>& vIn) {
  vector<double> conc;
  int size = vIn.size() * vIn[0].size();
  conc.reserve(size);
  for (int i = 0; i < vIn.size(); ++i)
    for (int j = 0; j < vIn[0].size(); ++j)
      conc[j*vIn.size() + i] = vIn[i][j];
  mat res(&conc.front(), vIn.size(), vIn[0].size(), false);
  return res;
}