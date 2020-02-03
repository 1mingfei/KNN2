#include "gbCnf.h"
#include "KNHome.h"
#include "FCCEmbededCluster.h"


void KNHome::createOrderedSingle(const int& i, \
                         gbCnf& cnfModifier, \
                         const vector<int>& dupFactors, \
                         const double& LC, \
                         const string& POT, \
                         const FCCEmbededCluster::occupInfo_256& o256) {

  Config c0 = cnfModifier.getFCCConv(LC, "Al", dupFactors);

  // cnfModifier.writeCfgData(c0, to_string(i) + ".cfg");

  vector<pair<string, string>> elemPairs = {
                                            {"Al", "Mg"}, \
                                            {"Al", "Zn"}, \
                                            {"Mg", "Zn"}, \
                                           };

  const vector<pair<int, int>>& jumpPairsRef = o256.jumpPairs[i];
  int index = 0, subIndex = 0;  // config and end index, respectively
  for (const auto& elemPair : elemPairs) {

    Config c1 = cnfModifier.embedCluster(c0, elemPair, o256, i);
    cnfModifier.writeCfgData(c1, "L10_" + to_string(index) + ".cfg");
    Config cS = c0;

    for (int j = 0; j < jumpPairsRef.size(); ++j) {

      string tmpType;
      if ( (j == 0) || \
           ((j > 0) && (jumpPairsRef[j].first != jumpPairsRef[j-1].first)) ) {
        // cout << j << " " << jumpPairsRef[j].first << endl;
        string baseDir = "config" + to_string(index) + "/s";
        string mkBaseDir = "mkdir -p " + baseDir;
        const char *cmkBaseDir = mkBaseDir.c_str();
        const int dir_err = std::system(cmkBaseDir);
        if (-1 == dir_err) {
          cout << "Error creating directory!\n";
          exit(1);
        }
        cS = c0;
        cS = cnfModifier.embedCluster(c1, elemPair, o256, i);

        tmpType = cS.atoms[jumpPairsRef[j].first].tp;
        cS.atoms[jumpPairsRef[j].first].tp = "X";

        cnfModifier.writeCfgData(cS, "config" + to_string(index) + \
                                     "/s/start.cfg");

        std::sort(cS.atoms.begin(), cS.atoms.end());
        map<string, int> elemName = cnfModifier.writePOSCAR(cS, \
                              "config" + to_string(index) + "/s/POSCAR");
        prepVASPFiles(baseDir, dupFactors, elemName, POT);
        subIndex = 0;
      }


      string subDir = "config" + to_string(index) \
                      + "/e_" + to_string(subIndex);
      const char *csubDir = subDir.c_str();
      const int dir_err = mkdir(csubDir, \
                                S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      if (-1 == dir_err) {
        cout << "Error creating directory!\n";
        exit(1);
      }

      Config cF = cnfModifier.swapPair(c0, jumpPairsRef[j]);
      cF = cnfModifier.embedCluster(cF, elemPair, o256, i);
      cF.atoms[jumpPairsRef[j].first].tp = "X";


      string name1 = "config" + to_string(index) \
                     + "/e_" + to_string(subIndex) + "/";
      cnfModifier.writeCfgData(cF, name1 + "end.cfg");

      std::sort(cF.atoms.begin(), cF.atoms.end());

      map<string, int> elemName = cnfModifier.writePOSCAR(cF, \
                                                  name1 + "POSCAR");
      cout << "config " << index << " end " << subIndex \
           << " pair: " << jumpPairsRef[j].first \
           << " " << jumpPairsRef[j].second << "\n";

      prepVASPFiles(name1, dupFactors, elemName, POT);

      ++subIndex;
      if ((j < jumpPairsRef.size()) \
          && (jumpPairsRef[j].first != jumpPairsRef[j+1].first)) {
        ++index;
      }
    }
  }
}

void KNHome::createOrdered(gbCnf& cnfModifier, \
                           const vector<int>& dupFactors,
                           const double& LC,
                           const string& POT) {

  FCCEmbededCluster::occupInfo_256 o256;

  // int quotient = NConfigs / nProcs;
  // int remainder = NConfigs % nProcs;
  // int nCycle = remainder ? (quotient + 1) : quotient;
  // for (int j = 0; j < nCycle; ++j) {
  //   for (int i = (j * nProcs); i < ((j + 1) * nProcs); ++i) {
  //     if ((i % nProcs != me) || (i >= NConfigs)) continue;
  //     createOrderedSingle(i, cnfModifier, dupFactors, LC, POT);
  //   }
  // }
  createOrderedSingle(0, cnfModifier, dupFactors, LC, POT, o256);

}

Config gbCnf::embedCluster(const Config& cIn, \
                           const pair<string, string>& Elems, \
                           const FCCEmbededCluster::occupInfo_256& o256, \
                           const int& i) {
  Config c0 = cIn;

  for (int j = 0; j < o256.mapping[i].size(); ++j) {
    if (o256.mapping[i][j] == 1) {
      c0.atoms[j].tp = Elems.first;
    } else if (o256.mapping[i][j] == 2) {
      c0.atoms[j].tp = Elems.second;
    }
  }
  // std::sort(c0.atoms.begin(), c0.atoms.end());

  return c0;
}