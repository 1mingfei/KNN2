#ifndef _GB_UTL_H_
#define _GB_UTL_H_

#include <numeric>
#include <string>
//#include <vector>

using std::cout;
using std::endl;
using std::string;
using std::vector;

inline vector<vector<double>> mtx33x33(const vector<vector<double>>& a,
                                       const vector<vector<double>>& b) {
  vector<vector<double>> c(3, vector<double>(3, 0));
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++) c[i][j] += a[i][k] * b[k][j];
  return c;
}

inline void printVec(const vector<double>& v) {
  for (int i = 0; i < 3; i++) cout << v[i] << " ";
  cout << endl;
}

inline double square11(double x) { return x * x; };
inline double square33(const vector<double>& v) {
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}
inline double innDot33(const vector<double>& a, const vector<double>& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
inline double relerr(double a, double b) { return fabs(a - b) / b; }
inline void scaleVec(vector<double>& v, double s) {
  for (auto &ee : v) ee *= s;
}
inline void crossProd33(const vector<double>& a, const vector<double>& b,
                        vector<double>& res) {
  res[0] = a[1] * b[2] - a[2] * b[1];
  res[1] = a[2] * b[0] - a[0] * b[2];
  res[2] = a[0] * b[1] - a[1] * b[0];
}
inline double vecInnProd33(const vector<double>& a, const vector<double>& b) {
  double res = 0;
  for (int i = 0; i < 3; i++) res += a[i] * b[i];
  return res;
}

// copy the string and return a char pointer
inline void split(const string& s, const char *delim, vector<string>& v) {
  char *dup = strdup(s.c_str());
  char *token = strtok(dup, delim);
  while (token != NULL) {
    v.push_back(string(token));
    // the call is treated as a subsequent calls to strtok:
    // the function continues from where it left in previous invocation
    token = strtok(NULL, delim);
  }
  free(dup);
}

inline void outMkdir(const string& mdir) {
  struct stat buf;
  if ((stat(mdir.c_str(), &buf) == 0) != 1) mkdir(mdir.c_str(), S_IRWXU);
}

#endif  // _GB_UTL_H_
