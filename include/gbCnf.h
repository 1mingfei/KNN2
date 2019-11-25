#ifndef _GB_CNF_H_
#define _GB_CNF_H_

#include "armadillo"
#include "KNHome.h"
#include <math.h>
using arma::mat;
using arma::vec;
using std::vector;
using std::pair;

class KNHome::gbCnf {
  KNHome& hm;
  unordered_map<string, string>& sparams;

public:
  double rcut;
  gbCnf(KNHome& x, double rc = 3.0)
      : hm(x),
        sparams(x.sparams), 
        rcut(rc) {};

  /* gbInCnf.cpp */
  Config readLmpData(const string&);
  Config readCfg(const string&);

  /* output
   * gbOutCnf.cpp */
  double findMass(string);
  void writeLmpData(Config&, string);
  //void writeLmpDataDebug(Config&, string);
  void writeCfgData(const Config& c, string);
  void writePOSCARVis(Config&, string, string);
  void writePOSCAR(Config&, string);

  /* gbInCnf.cpp */
  void cnvVec2Mat(const vector<double>&, Config& c);
  void cnvMat2Vec(Config&);
  vector<double> cnvVecXY2VecAng(const vector<double>& v);

  /* neighbor.cpp*/
  int getExpdParam(const Config&, const double);
  vector<KNAtom> expandCellZ(const Config&, const int);
  double calDist(const vector<double>, const KNAtom&, const KNAtom&);
  double calDistPrl(const vector<double>&, const KNAtom&, const KNAtom&);
  void getNBL(Config&, double);

  /* gbBox.cpp */
  void initBox(Config&);
  void wrapAtomPos(Config&);
  void wrapAtomPrl(Config&);
  void cnvprl2pst(Config&);
  void cnvpst2prl(Config&);

  /* FCCConfig.cpp */
  Config getFCCConv(const double, const string, const vector<int>&);

  /* KNSolidSol.cpp
   * elems stores name of elements
   * nums stores their corresponding numbers
   */
  void getRandConf(Config& cnf, const vector<string>& elems,\
                   const vector<int>& nums);

  void getRandConfUniformDist(Config& cnf, vector<string>& elems,\
                              const vector<int>& nums);
  vector<pair<int, int>> getPairToSwap(Config&);
  Config swapPair(const Config&, pair<int, int>);

  /* KNEncodeCart.cpp */
  vector<vector<string>> encodeConfig(Config&, const vector<int>&, const double,\
                                      vector<string>&, const vector<int>&,\
                                      const bool);

  vector<pair<int, int>> readPairs(string);
  Config rotateJumpPair(Config&, const vector<int>, const Config&);
  Config rotateConfig(Config& cfgOld, const vector<double>& v2);

  vec getCenterShift(Config&);
  void shiftToCenter(Config&, vector<double>&);
  mat getJumpCoor(const Config&, const vector<int>, const Config&);
                       
  /* KNBondCount.cpp */
  map<string, int> countPairs(Config&, const vector<string>&, \
                                       const vector<int>& );

};

#include "Elem.inl"
#endif
