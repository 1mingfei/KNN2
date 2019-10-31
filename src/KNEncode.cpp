/*
 * Author: 1mingfei 
 * Date:   2019-05-25
 * Purpose: encoding neighbor atoms of the jumping pairs
 * self-explained
 */

#include "gbCnf.h"
#include "KNHome.h"
#include "armadillo"
using arma::mat;
using arma::vec;


/*
 * Sort points lexicographically 
 * this will only work for the on-lattice cases
 * because it has no buffer range when comparing
 */
inline void sortAtomLexi(vector<KNAtom>& atmList) {
  sort(
    atmList.begin(), atmList.end(),
    [](const KNAtom& a, const KNAtom& b) -> bool
      { return (a.pst[X] < b.pst[X]) || 
               (a.pst[X] == b.pst[X] && a.pst[Y] < b.pst[Y]) ||
               (a.pst[X] == b.pst[X] && a.pst[Y] == b.pst[Y] && 
                a.pst[Z] < b.pst[Z]); }
    );
}

/*
 * Sort input vectors by their first two element
 * e.g. 
 * before sorting:       after sorting:
 * 1 2 4 9               0 2 8 9
 * 2 9 7 3         ==>   0 3 3 4
 * 0 3 3 4         ==>   1 2 4 9
 * 0 2 8 9               2 9 7 3
 */
inline bool sortPairs(const vector<int>& lhs, const vector<int>& rhs) {
  assert(lhs.size() >= 2 && rhs.size() >= 2);
  if (lhs[0] == rhs[0]) {
    return lhs[1] < rhs[1];
  } else {
    return lhs[0] < rhs[0];
  }
}

/* output vectors to file */
template<class T>
inline void writeVector(const string& fname, \
                        const vector<T>& v) {
  ofstream ofs(fname, std::ofstream::app);
  for (const auto& val : v) {
    ofs << std::setw(4) << val << " ";
  }
  ofs << endl;
}

/* output vectors to file */
template<class T>
inline void writeVector(const string& fname, const int i, const int j, \
                        const vector<T>& v) {
  ofstream ofs(fname, std::ofstream::app);
  ofs << "config " << i << " end " << j << " ";
  for (const auto& val : v) {
    ofs << std::setw(4) << val << " ";
  }
  ofs << endl;
}

inline mat calculateRotateMatrix(
    const vec& A, const vec& B) {
  mat R; R.eye(3, 3);
  vec v = cross(A, B);
  double c = dot(A, B);
  mat vx(3, 3);
  vx << 0.0 << -v[2] << v[1] << arma::endr 
     << v[2] << 0.0 << -v[0] << arma::endr
     << -v[1] << v[0] << 0.0 <<arma::endr;

  R = R + vx + vx * vx / (1.0 + c);
  return R;
}

inline vector<double> rotateMatrix(
    const mat& A, const vector<double>& X) {
  vec Xv(3);
  Xv << X[0] << X[1] << X[2];
  vec B = A * Xv;
  return {B[0], B[1], B[2]};
}

inline vector<double> rotateMatrix(
    const mat& A, const vector<double>& X, const vec& center) {
  vec Xv(3);
  Xv << X[0] - center[0] << X[1] - center[1] << X[2] - center[2];
  vec B = A * Xv;
  //return {B[0] + center[0], B[1] + center[1], B[2] + center[2]};
  return {B[0], B[1], B[2]};
}

vec KNHome::gbCnf::getCenterShift(Config& c) {
  //unwrap atoms
  vector<double> ave(3, 0.0);
  wrapAtomPrl(c);
  for (auto&& atm : c.atoms) {
    for (int i = 0; i < 3; ++i) {
      if (atm.prl[i] >= 0.5)
        atm.prl[i] -= 1.0;
      ave[i] += atm.prl[i];
    }
  }
  for (int i = 0; i < 3; ++i) {
    ave[i] /= c.atoms.size();
  }
  vec res(3);
  res << ave[0] - 0.5 << ave[1] - 0.5 << ave[2] - 0.5;
  return res;
}

void KNHome::KNEncode() {
  gbCnf cnfModifier(*this);
  vector<string> elems = vsparams["elems"];
  int NConfigs = iparams["NConfigs"];
  int NBars = iparams["NBarriers"];
  double RCut = dparams["RCut"];
  //string fpname = "pairs.txt";//
  string fpname = sparams["PairFile"];
  vector<vector<int>> pairs = readPairs(fpname);
  MPI_Barrier(MPI_COMM_WORLD);
  if (me == 0)
    std::cout << "encoding configurations\n";
  //int quotient = NConfigs / nProcs;
  //int remainder = NConfigs % nProcs;
  int quotient = pairs.size() / nProcs;
  int remainder = pairs.size() % nProcs;
  int nCycle = remainder ? (quotient + 1) : quotient;
  for (int j = 0; j < nCycle; ++j) {
    for (int i = (j * nProcs); i < ((j + 1) * nProcs); ++i) {
      if ((i % nProcs != me) || (i >= pairs.size())) continue;
      //string fname = "in.cfg";//
      string fname = "config" + to_string(pairs[i][0]) + "/s/start.cfg";

      Config cfg = cnfModifier.readCfg(fname);
      vector<int> pair = {pairs[i][2], pairs[i][3]};
      vector<string> codes;
      
      vector<int> resId = cnfModifier.encodeConfig(cfg, pair, RCut, codes);

      //writeVector<int>(to_string(i) + ".txt", pairs[i][0], pairs[i][1],\
      //                 resId);
      writeVector<string>("encode.out.txt", pairs[i][0], pairs[i][1],\
                         codes);
      /*
      vector<int> tmpId = cnfModifier.encodeConfig(cfg, pair, RCut, codes);
      writeVector<int>(to_string(i) + ".txt", tmpId);
      writeVector<string>(to_string(i) + ".txt", codes);

      Config cNew = cfg;
      cNew.atoms.clear();
      cNew.natoms = tmpId.size();
      for (unsigned int i = 0 ; i < tmpId.size(); ++i) {
        cNew.atoms.push_back(cfg.atoms[tmpId[i]]);
      }
#ifdef DEBUG
      cnfModifier.writeCfgData(cNew, to_string(i) + "_encode_before.cfg");
#endif

      Config cfgRotated = cnfModifier.rotate(cNew, pair);
      vector<string> codesR;
      vector<int> resId;
      sortAtomLexi(cfgRotated.atoms);
      for (const auto& atm : cfgRotated.atoms) {
        resId.push_back(atm.id);
        codesR.push_back(atm.tp);
      }
      //vector<int> resId = cnfModifier.encodeConfig(cfgRotated, pair, RCut,\
      //                                             codesR);
      writeVector<int>(to_string(i) + ".txt", resId);
      writeVector<string>(to_string(i) + ".txt", codesR);

#ifdef DEBUG
      cnfModifier.writeCfgData(cfgRotated, to_string(i) + "_encode_after.cfg");
#endif
      */

    }
  }
}


vector<vector<int>> KNHome::readPairs(const string& fname) {
  ifstream ifs(fname.empty() ? sparams["datafile"] : fname, std::ifstream::in);
  vector<vector<int>> res;
  vector<string> s;
  string buff;
  getline(ifs, buff); //for comment line Config Barrier atom1 atom2
  while (getline(ifs, buff)) {
    s.clear();
    split(buff, " ", s);
    if (s[0] == "Delta") //for the very last line
      continue;
    vector<int> tmp;
    //assert(s.size() == 4);
    tmp.push_back(stoi(s[1]));
    tmp.push_back(stoi(s[3]));
    for (int i = 5; i < 7; ++i) {
      tmp.push_back(stoi(s[i]));
    }
    res.push_back(tmp);
  }
  sort(res.begin(), res.end(), sortPairs);
  return res;
}

vector<int> KNHome::gbCnf::encodeConfig(Config& cnf,
                                        const vector<int> pair, \
                                        double RCut, vector<string>& codes) {
  assert(pair.size() == 2); //the size of input pair must equals 2
  getNBL(cnf, RCut);

  /* reverse wrap back so the lexi order is always not affected by Periodic 
   * boundary conditions */
  vector<KNAtom> atmList;
  vector<int> tmpId;
  set<int> hsh;
  for (int i = 0; i < 2; ++i) {
    KNAtom atm = cnf.atoms[pair[i]];
    for (const int ii : atm.NBL) {
      KNAtom nbAtm = cnf.atoms[ii];
      for (int j = 0; j < DIM; ++j) {
        if (std::abs(atm.prl[j] - 1.0) < RCut / cnf.length[j]) {
          if ((1.0 - nbAtm.prl[j]) < RCut / cnf.length[j]) {
            nbAtm.prl[j] -= 1.0;
          } else if ((nbAtm.prl[j] - 1.0) < RCut / cnf.length[j]) {
            nbAtm.prl[j] += 1.0;
          }
        }
      }
      if (hsh.find(nbAtm.id) == hsh.end()) {
        hsh.insert(nbAtm.id);
        if ((nbAtm.id != cnf.atoms[pair[0]].id) && \
            (nbAtm.id != cnf.atoms[pair[1]].id)) {
          atmList.push_back(nbAtm);
        }
      }
    }
  }
  //sortAtomLexi(atmList);
  vector<string> tmpCodes;
  for (const auto& atm : atmList) {
    tmpId.push_back(atm.id);
    tmpCodes.push_back(atm.tp);
  }

#ifdef DEBUG
  writeVector<int>("debug_encode.txt", tmpId);
  writeVector<string>("debug_encode.txt", tmpCodes);
#endif

  Config cNew = cnf;
  cNew.atoms.clear();
  cNew.natoms = tmpId.size();
  for (unsigned int i = 0 ; i < tmpId.size(); ++i) {
    cNew.atoms.push_back(cnf.atoms[tmpId[i]]);
  }
#ifdef DEBUG
  writeCfgData(cNew, "debug_encode_before.cfg");
#endif

#ifdef DEBUG
  cout << "before Rotate\n";
  for (int i = 0; i < 5; ++i)
    cout << cNew.atoms[i].prl[0] << " " << cNew.atoms[i].prl[1] << " " << cNew.atoms[i].prl[2] << "\n";
  cout << "\n";
#endif 

  Config cfgRotated = rotate(cNew, pair);
  vector<int> resId;
  codes.push_back(cnf.atoms[pair[1]].tp);//
  sortAtomLexi(cfgRotated.atoms);
  for (const auto& atm : cfgRotated.atoms) {
    resId.push_back(atm.id);
    codes.push_back(atm.tp);
  }
  //vector<int> resId = cnfModifier.encodeConfig(cfgRotated, pair, RCut,\
  //                                             codesR);
#ifdef DEBUG
  writeCfgData(cfgRotated, "debug_encode_after.cfg");
#endif
  return resId;
}

Config KNHome::gbCnf::rotate(Config& cnf, const vector<int> pair) {

#ifdef DEBUG
  cout << "before rotate\n";
  for (int i = 0; i < 5; ++i)
    cout << cnf.atoms[i].prl[0] << " " << cnf.atoms[i].prl[1] << " " << cnf.atoms[i].prl[2] << "\n";
  cout << "\n";
#endif 

  cnvprl2pst(cnf);
  wrapAtomPos(cnf);
#ifdef DEBUG
  writeCfgData(cnf, "debug_1.cfg");
#endif
  cnvpst2prl(cnf);
#ifdef DEBUG
  writeCfgData(cnf, "debug_2.cfg");
#endif

// #ifdef DEBUG
//   cout << "rotate step3\n";
//   for (int i = 0; i < 5; ++i)
//     cout << cnf.atoms[i].prl[0] << " " << cnf.atoms[i].prl[1] << " " << cnf.atoms[i].prl[2] << "\n";
//   cout << "\n";
// #endif 

  int id1 = pair[0], id2 = pair[1];
  vector<double> v1(3, 0.0);
  for (unsigned int i = 0; i < DIM; ++i) {
    v1[i] = cnf.atoms[id2].prl[i] - cnf.atoms[id1].prl[i];
  }
  vec v1a(3);
  v1a << v1[X] << v1[Y] << v1[Z];
  v1a = arma::normalise(v1a);
  vector<double> v2 = {1.0, 1.0, 0.0};
  vec v2a;
  v2a << v2[X] << v2[Y] << v2[Z];
  v2a = arma::normalise(v2a);
  mat R = calculateRotateMatrix(v1a, v2a);
  cout << "initial jump direction: \n";
  v1a.print();
  cout << "reference jump direction: \n";
  v2a.print();
  cout << "rotation matrix: \n";
  R.print();
  Config res = cnf;

  //the distance of center of atoms to the center of the cell
  vec centerShift = getCenterShift(cnf); 
  centerShift.print();

  for (auto&& atm : res.atoms) {
    vector<double> prl({atm.prl[0], atm.prl[1], atm.prl[2]});
    prl = rotateMatrix(R, prl, centerShift);
    atm.prl[0] = prl[0], atm.prl[1] = prl[1], atm.prl[2] = prl[2];
  }
#ifdef DEBUG
  writeCfgData(res, "debug_3.cfg");
#endif
  //cnvprl2pst(res);
  wrapAtomPrl(res);
#ifdef DEBUG
  writeCfgData(res, "debug_4.cfg");
#endif
  //cnvpst2prl(res);
  return res;
}
