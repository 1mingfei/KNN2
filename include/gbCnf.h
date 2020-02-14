#ifndef _GB_CNF_H_
#define _GB_CNF_H_

#include "armadillo"
#include "KNHome.h"
#include "FCCEmbededCluster.h"
#include <math.h>

using arma::mat;
using arma::vec;
using std::vector;
using std::pair;
using std::map;
using std::unordered_multimap;

using keras2cpp::Model;
using keras2cpp::Tensor;

class gbCnf {
  // KNHome& hm;
  unordered_map<string, string>& sparams;

public:
  double rcut;
  gbCnf(unordered_map<string, string>&);

  /* gbInCnf.cpp */
  Config readLmpData(const string&);
  Config readCfg(const string&);
  Config readPOSCAR(const string&);

  /* output
   * gbOutCnf.cpp */
  double findMass(string);
  void writeLmpData(Config&, string);
  //void writeLmpDataDebug(Config&, string);
  void writeCfgData(const Config&, string);
  void writeCfgAux(const Config&, const vector<int>&, string);

  void writePOSCARVis(Config&, string, string);
  map<string, int> writePOSCAR(Config&, string);

  /* gbInCnf.cpp */
  void cnvVec2Mat(const vector<double>&, Config&);
  void cnvMat2Vec(Config&);
  vector<double> cnvVecXY2VecAng(const vector<double>&);

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
  void perturb(Config&);

  /* FCCConfig.cpp */
  Config getFCCConv(const double, const string, const vector<int>&);

  /* KNSolidSol.cpp
   * elems stores name of elements
   * nums stores their corresponding numbers
   */
  void getRandConf(Config&, \
                   std::default_random_engine&, \
                   const vector<string>&, \
                   const vector<int>&);

  void getRandConfUniformDist(Config&, \
                              std::default_random_engine&, \
                              vector<string>&, \
                              const vector<int>&);

  vector<pair<int, int>> getPairToSwap(Config&);

  Config swapPair(const Config&, pair<int, int>);

  /* KNEncodeCart.cpp */
  vector<vector<string>> encodeConfig(Config&, \
                                      const vector<int>&, \
                                      const double, \
                                      vector<string>&, \
                                      const vector<int>&, \
                                      const bool);

  vector<pair<int, int>> readPairs(string);
  Config rotateJumpPair(Config&, const vector<int>, const Config&);
  Config rotateConfig(Config&, const vector<double>&);

  vec getCenterShift(Config&);
  void shiftToCenter(Config&, vector<double>&);
  mat getJumpCoor(const Config&, const vector<int>, const Config&);

  /* KNBondCount.cpp */
  map<string, int> countPairs(Config&, const vector<string>&, \
                                       const vector<int>& );

  /* KMCSimulation.cpp */
  vector<double> calBarrierAndEdiff(Config&, \
                                   const double&, \
                                   const double&, \
                                   const string&, \
                                   unordered_map<string, double>&, \
                                   Model&, \
                                   Model&, \
                                   const pair<int, int>&);

  /* findClusters.cpp */

  // This function give a set which contains all the solute atoms ID
  unordered_set<int> findSoluteAtoms(const Config&, const string&);

  // This is a helper function which is to find X clusters after
  // removing a certain element and returns X.
  int helperBFS(const Config&, \
                const unordered_set<int>&, \
                unordered_multimap<int, int>&, \
                map<int, int>&);

  void getLargestClts(const int&, \
                      const int&, \
                      unordered_multimap<int, int>&, \
                      map<int, int>&);

  void helperAddFNNs(const Config&, \
                     unordered_multimap<int, int>&, \
                     map<int, int>&);
  // This function returns X largest clusters with FNBs
  map<int, int> findAtm2Clts(Config&, const int&, const string&);

  /* KNOrdered.cpp */
  Config embedCluster(const Config&, \
                      const pair<string, string>&, \
                      const FCCEmbededCluster::AFOccupInfo_256&, \
                      const int&);
};

#include "Elem.inl"
#endif
