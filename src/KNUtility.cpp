#include "KNUtility.h"

mat vvd2mat(vvd& mIn) {
  vd conc;
  int size = mIn.size() * mIn[0].size();
  conc.reserve(size);
  for (int i = 0; i < mIn.size(); ++i)
    for (int j = 0; j < mIn[0].size(); ++j)
      conc[j*mIn.size() + i] = mIn[i][j];
  mat res(&conc.front(), mIn.size(), mIn[0].size(), false);
  return res;
}

vec vd2vec(vd& vIn) {
  vec res(&vIn.front(), vIn.size(), false);
  return res;
}

vd mat2vd(mat& mIn) {
  vd res(mIn.n_rows * mIn.n_cols, 0.0);
  for (int i = 0; i < mIn.n_rows; ++i)
    for (int j = 0; j < mIn.n_cols; ++j)
      res[i * mIn.n_cols + j] = mIn.at(i, j);
  return res;
}