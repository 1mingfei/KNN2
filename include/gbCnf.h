#ifndef _GB_CNF_H
#define _GB_CNF_H

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
  double calDistPrl(const vector<double>, const KNAtom&, const KNAtom&);
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

  /* KNEncode.cpp */
  vector<int> encodeConfig(Config&, const vector<int>, double, vector<string>&);
  vector<pair<int, int>> readPairs(string);
  Config rotate(Config&, const vector<int>, const Config&);
  vec getCenterShift(Config&);
  void shiftToCenter(Config&, vector<double>&);
  /* KNBondCount.cpp */
  map<string, int> countPairs(Config&, const vector<string>&, \
                                       const vector<int>& );

};

#include "Elem.inl"
#endif
