/*
 * @Author: chaomy
 * @Date:   2018-06-20 15:52:26
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-11-11 13:51:06
 */

#include "gbCnf.h"

/**************************************************
 * convert vector cell to matrix cell
 **************************************************/
void KNHome::gbCnf::cnvVec2Mat(const vector<double>& v, Config& c) {
  c.bvx[0] = v[0];
  c.bvy[0] = v[1] * cos(v[5]);
  c.bvy[1] = v[1] * sin(v[5]);
  c.bvz[0] = v[2] * cos(v[4]);
  c.bvz[1] = v[2] * cos(v[3]) * sin(v[5]) -
             ((v[2] * cos(v[4]) - v[2] * cos(v[3]) * cos(v[5])) / tan(v[5]));
  c.bvz[2] = sqrt(v[2] * v[2] - c.bvz[0] * c.bvz[0] - c.bvz[1] * c.bvz[1]);
}

/**************************************************
 * convert matrix cell to vector cell
 **************************************************/
void KNHome::gbCnf::cnvMat2Vec(Config& c) {
  vector<double>& v = c.cell = vector<double>(9, 0);
  v[3] = std::sqrt(innDot33(c.bvx, c.bvx));
  v[4] = std::sqrt(innDot33(c.bvy, c.bvy));
  v[5] = std::sqrt(innDot33(c.bvz, c.bvz));
  v[6] = std::acos(innDot33(c.bvx, c.bvy) / (v[X + 3] * v[Y + 3]));
  v[7] = std::acos(innDot33(c.bvx, c.bvz) / (v[X + 3] * v[Z + 3]));
  v[8] = std::acos(innDot33(c.bvy, c.bvz) / (v[Y + 3] * v[Z + 3]));
}

/**************************************************
 * Convert vector xyz to angle alpha, beta gamma
 * lx ly lz xy xz yz
 * 0  1  2  3  4  5
 **************************************************/
vector<double> KNHome::gbCnf::cnvVecXY2VecAng(const vector<double>& v) {
  // a b c alpha beta gamma
  // 0 1 2 3     4    5
  vector<double> r(6, 0);
  r[0] = v[0];
  r[1] = std::sqrt(square11(v[1]) + square11(v[3]));
  r[2] = std::sqrt(square11(v[2]) + square11(v[4]) + square11(v[5]));
  r[3] = std::acos((v[3] * v[4] + v[1] * v[5]) / (r[1] * r[2]));
  r[4] = std::acos(v[4] / r[2]);
  r[5] = std::acos(v[3] / r[1]);
  return r;
}

/**************************************************
 * read lammps data files
 **************************************************/
/*
Config KNHome::gbCnf::readLmpData(const string& fname) {
  ifstream ifs(fname.empty() ? sparams["datafile"] : fname, std::ifstream::in);
  string buff;
  Config cnf;
  vector<string> s;
  unordered_map<string, int> mp;
  getline(ifs, buff);
  sscanf(buff.c_str(), "# %lf %lf", &cnf.oldEngy, &cnf.engy);

  getline(ifs, buff, ' ');
  cnf.natoms = stoi(buff);
  getline(ifs, buff);

  getline(ifs, buff, ' ');  // atom types
  cnf.ntypes = stoi(buff);
  getline(ifs, buff);

  vector<double> v(6);
  for (auto&& i : {0, 1, 2}) {  // lo hi
    getline(ifs, buff);
    sscanf(buff.c_str(), "%lf %lf", &cnf.cell[i], &cnf.cell[i + 3]);
    v[i] = cnf.length[i] = cnf.cell[3 + i] - cnf.cell[i];
    cnf.center[i] = 0.5 * (cnf.cell[3 + i] + cnf.cell[i]);
  }
  getline(ifs, buff);
  s.clear();
  split(buff, " ", s);
  cnf.cell[6] = stof(s[0]), cnf.cell[7] = stof(s[1]), cnf.cell[8] = stof(s[2]);
  cnvVec2Mat(cnvVecXY2VecAng(v), cnf);

  while (getline(ifs, buff)) {
    s.clear();
    split(buff, " ", s);
    if (s.size() > 3) break;
  }
  for (int i = 0; i < cnf.natoms; ++i) {
    Atom a;
    sscanf(buff.c_str(), "%d %d %lf %lf %lf", &a.id, &a.tp, &a.pst[0],
           &a.pst[1], &a.pst[2]);
    a.id -= 1;  // change from 1 based to 0 based
    cnf.atoms.push_back(a);
    getline(ifs, buff);
  }
  sort(cnf.atoms.begin(), cnf.atoms.end());  // sort based on atom id
  return cnf;
}
*/

