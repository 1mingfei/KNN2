#include "gbCnf.h"
inline vector<string> splitStr2Str(const string& stringIn) {
  // Used to split string around spaces.
  stringstream ss(stringIn);
  vector<string> wordVector;
  // Traverse through all words
  do {
      // Read a word
      string word;
      ss >> word;
      wordVector.push_back(word);

  // While there is more to read
  } while (ss);
  return wordVector;
}

inline vector<int> splitStr2Int(const string& stringIn) {
  // Used to split string around spaces.
  stringstream ss(stringIn);
  vector<int> wordVector;
  // Traverse through all words
  do {
      // Read a word
      string word;
      ss >> word;
      wordVector.push_back(stoi(word));

  // While there is more to read
  } while (ss);
  return wordVector;
}

gbCnf::gbCnf(int& meIn, int& nProcsIn)
      : me(meIn), \
        nProcs(nProcsIn), \
        rcut(3.0) {};

/**************************************************
 * convert vector cell to matrix cell
 **************************************************/
void gbCnf::cnvVec2Mat(const vector<double>& v, Config& c) {
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
void gbCnf::cnvMat2Vec(Config& c) {
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
vector<double> gbCnf::cnvVecXY2VecAng(const vector<double>& v) {
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
 * read cfg data files
 **************************************************/

Config gbCnf::readCfg(const string& fname) {
  ifstream ifs(fname, std::ifstream::in);
  string buff;
  Config cnf;
  getline(ifs, buff);
  sscanf(buff.c_str(), "Number of particles = %i", &cnf.natoms);
  getline(ifs, buff); // A = 1.0 Angstrom (basic length-scale)
  getline(ifs, buff);
  sscanf(buff.c_str(), "H0(1,1) = %lf A", &cnf.bvx[X]);
  getline(ifs, buff);
  sscanf(buff.c_str(), "H0(1,2) = %lf A", &cnf.bvx[Y]);
  getline(ifs, buff);
  sscanf(buff.c_str(), "H0(1,3) = %lf A", &cnf.bvx[Z]);
  getline(ifs, buff);
  sscanf(buff.c_str(), "H0(2,1) = %lf A", &cnf.bvy[X]);
  getline(ifs, buff);
  sscanf(buff.c_str(), "H0(2,2) = %lf A", &cnf.bvy[Y]);
  getline(ifs, buff);
  sscanf(buff.c_str(), "H0(2,3) = %lf A", &cnf.bvy[Z]);
  getline(ifs, buff);
  sscanf(buff.c_str(), "H0(3,1) = %lf A", &cnf.bvz[X]);
  getline(ifs, buff);
  sscanf(buff.c_str(), "H0(3,2) = %lf A", &cnf.bvz[Y]);
  getline(ifs, buff);
  sscanf(buff.c_str(), "H0(3,3) = %lf A", &cnf.bvz[Z]);
  getline(ifs, buff); // .NO_VELOCITY.
  getline(ifs, buff);
  int entry = 3;
  sscanf(buff.c_str(), "entry_count = %i", &entry);
  vector<string> s;

  for (int i = 0; i < cnf.natoms; ++i) {
    KNAtom a;
    double mass;
    getline(ifs, buff);
    sscanf(buff.c_str(), "%lf", &mass);
    getline(ifs, buff);
    s.clear();
    split(buff, " ", s);
    a.tp = s[0];
    getline(ifs, buff);
    sscanf(buff.c_str(), "%lf %lf %lf", &a.prl[0], &a.prl[1], &a.prl[2]);
    a.id = i;
    cnf.atoms.push_back(a);
  }

  cnf.length[X] = cnf.bvx[X];
  cnf.length[Y] = cnf.bvy[Y];
  cnf.length[Z] = cnf.bvz[Z];

  // std::sort(cnf.atoms.begin(), cnf.atoms.end());

  return cnf;
}

Config gbCnf::readPOSCAR(const string& fname) {
  ifstream ifs(fname, std::ifstream::in);
  string buff;
  Config cnf;
  getline(ifs, buff); // comment line
  getline(ifs, buff); // scale factor, I skip this since I will always use 1
  getline(ifs, buff);
  sscanf(buff.c_str(), "%lf %lf %lf", &cnf.bvx[X], &cnf.bvx[Y], &cnf.bvx[Z]);
  getline(ifs, buff);
  sscanf(buff.c_str(), "%lf %lf %lf", &cnf.bvy[X], &cnf.bvy[Y], &cnf.bvy[Z]);
  getline(ifs, buff);
  sscanf(buff.c_str(), "%lf %lf %lf", &cnf.bvz[X], &cnf.bvz[Y], &cnf.bvz[Z]);
  getline(ifs, buff);

  string elem_names_string;

  vector<string> elem_names;
  split(buff, " ", elem_names);

  getline(ifs, buff);

  vector<string> elem_counts_str;
  split(buff, " ", elem_counts_str);
  vector<int> elem_counts;
  for (int i = 0; i < elem_counts_str.size(); ++i) {
    elem_counts.push_back(stoi(elem_counts_str[i]));
  }

  int sum = accumulate(elem_counts.begin(), elem_counts.end(), 0);

  cnf.natoms = sum;

  getline(ifs, buff);
  if (buff[0] == 'D' || buff[0] == 'd') {
    int count = 0;
    for (int i = 0; i < elem_counts.size(); ++i) {
      for (int j = 0; j < elem_counts[i]; ++j) {
        KNAtom a;
        getline(ifs, buff);
        sscanf(buff.c_str(), "%lf %lf %lf", &a.prl[0], &a.prl[1], &a.prl[2]);
        a.id = count++;
        a.tp = elem_names[i];
        cnf.atoms.push_back(a);
      }
    }
    cnvprl2pst(cnf);
  } else {
    int count = 0;
    for (int i = 0; i < elem_counts.size(); ++i) {
      for (int j = 0; j < elem_counts[i]; ++j) {
        KNAtom a;
        getline(ifs, buff);
        sscanf(buff.c_str(), "%lf %lf %lf", &a.pst[0], &a.pst[1], &a.pst[2]);
        a.id = count++;
        a.tp = elem_names[i];
        cnf.atoms.push_back(a);
      }
    }
    cnvpst2prl(cnf);
  }
  cnf.length[X] = cnf.bvx[X];
  cnf.length[Y] = cnf.bvy[Y];
  cnf.length[Z] = cnf.bvz[Z];
  return cnf;
}

/**************************************************
 * read cfg data files with auxiliary
 **************************************************/

Config gbCnf::readCfgCluster(const string& fname, \
                             vector<unordered_set<int>>& map) {
  ifstream ifs(fname, std::ifstream::in);
  string buff;
  Config cnf;
  getline(ifs, buff);
  sscanf(buff.c_str(), "Number of particles = %i", &cnf.natoms);
  getline(ifs, buff); // A = 1.0 Angstrom (basic length-scale)
  getline(ifs, buff);
  sscanf(buff.c_str(), "H0(1,1) = %lf A", &cnf.bvx[X]);
  getline(ifs, buff);
  sscanf(buff.c_str(), "H0(1,2) = %lf A", &cnf.bvx[Y]);
  getline(ifs, buff);
  sscanf(buff.c_str(), "H0(1,3) = %lf A", &cnf.bvx[Z]);
  getline(ifs, buff);
  sscanf(buff.c_str(), "H0(2,1) = %lf A", &cnf.bvy[X]);
  getline(ifs, buff);
  sscanf(buff.c_str(), "H0(2,2) = %lf A", &cnf.bvy[Y]);
  getline(ifs, buff);
  sscanf(buff.c_str(), "H0(2,3) = %lf A", &cnf.bvy[Z]);
  getline(ifs, buff);
  sscanf(buff.c_str(), "H0(3,1) = %lf A", &cnf.bvz[X]);
  getline(ifs, buff);
  sscanf(buff.c_str(), "H0(3,2) = %lf A", &cnf.bvz[Y]);
  getline(ifs, buff);
  sscanf(buff.c_str(), "H0(3,3) = %lf A", &cnf.bvz[Z]);
  getline(ifs, buff); // .NO_VELOCITY.
  getline(ifs, buff);
  int entry = 3;
  sscanf(buff.c_str(), "entry_count = %i", &entry);
  for (int i = 0; i < entry - 3; ++i) // skipping line name
    getline(ifs, buff);

  vector<string> s;

  int maxClusterID = -1;
  for (int i = 0; i < cnf.natoms; ++i) {
    KNAtom a;
    double mass;
    getline(ifs, buff);
    sscanf(buff.c_str(), "%lf", &mass);
    getline(ifs, buff);
    s.clear();
    split(buff, " ", s);
    a.tp = s[0];
    getline(ifs, buff);
    sscanf(buff.c_str(), "%lf %lf %lf %i", \
           &a.prl[0], &a.prl[1], &a.prl[2], &a.clusterID);
    a.id = i;
    maxClusterID = std::max(maxClusterID, a.clusterID);
    cnf.atoms.push_back(a);
  }

  cnf.length[X] = cnf.bvx[X];
  cnf.length[Y] = cnf.bvy[Y];
  cnf.length[Z] = cnf.bvz[Z];

  map.resize(maxClusterID + 1);

  for (int i = 0; i < cnf.natoms; ++i) {
    int clusterID = cnf.atoms[i].clusterID;
    map[clusterID].insert(i);
  }

  return cnf;
}