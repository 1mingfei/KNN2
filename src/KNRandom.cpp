#include "gbCnf.h"
#include "KNHome.h"
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

void KNHome::createRandom(gbCnf& cnfModifier, \
                          const int& NConfigs, \
                          const int& NBars,  \
                          const string& POT, \
                          const double& LC, \
                          const vector<int>& dupFactors, \
                          const vector<string>& elems, \
                          const vector<int>& nums) {

  int quotient = NConfigs / nProcs;
  int remainder = NConfigs % nProcs;
  int nCycle = remainder ? (quotient + 1) : quotient;
  for (int j = 0; j < nCycle; ++j) {
    for (int i = (j * nProcs); i < ((j + 1) * nProcs); ++i) {
      if ((i % nProcs != me) || (i >= NConfigs)) continue;

      Config c0 = cnfModifier.getFCCConv(LC, elems[0], dupFactors);
      std::default_random_engine rng(std::random_device{}() + me);

      cnfModifier.getRandConf(c0, rng, elems, nums);
      Config c0copy = c0;

      string baseDir = "config" + to_string(i) + "/s";
      string mkBaseDir = "mkdir -p " + baseDir;
      const char *cmkBaseDir = mkBaseDir.c_str();
      const int dir_err = std::system(cmkBaseDir);
      if (-1 == dir_err) {
        cout << "Error creating directory!\n";
        exit(1);
      }

      cnfModifier.writeCfgData(c0, "config" + to_string(i) + \
                                   "/s/start.cfg");

      // add small perturbation to break perfect fcc symmetry
      // this method is about to increase the chance
      // to find lower ground states for VASP software
      cnfModifier.cnvprl2pst(c0);
      cnfModifier.perturb(c0);

      map<string, int> elemName = cnfModifier.writePOSCAR(c0, \
                            "config" + to_string(i) + "/s/POSCAR");
      prepVASPFiles(baseDir, dupFactors, elemName, POT);
      vector<pair<int, int>> pairs = cnfModifier.getPairToSwap(c0copy);

      std::shuffle(std::begin(pairs), std::end(pairs), rng);

      int end = MIN(NBars, pairs.size());
      for (unsigned int k = 0; k < end; ++k) {

        string subDir = "config" + to_string(i) + "/e_" + to_string(k);
        const char *csubDir = subDir.c_str();
        const int dir_err = mkdir(csubDir, \
                                  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (-1 == dir_err) {
          cout << "Error creating directory!\n";
          exit(1);
        }

        Config c1 = cnfModifier.swapPair(c0copy, pairs[k]);
        string name1 = "config" + to_string(i) + "/e_" + to_string(k) + "/";
        cnfModifier.writeCfgData(c1, name1 + "end.cfg");

        // add small perturbation to break perfect fcc symmetry
        // this method is about to increase the chance
        // to find lower ground states for VASP software
        cnfModifier.cnvprl2pst(c1);
        cnfModifier.perturb(c1);

        map<string, int> elemName = cnfModifier.writePOSCAR(c1, \
                                                    name1 + "POSCAR");

        ofstream ofs("log.txt", std::ofstream::app);

        ofs << "config " << i << " end " << k \
            << " pair: " << pairs[k].first \
            << " " << pairs[k].second << "\n";

        ofs.close();

        prepVASPFiles(name1, dupFactors, elemName, POT);
      }
    }
  }
}

void KNHome::createRandomUniform(gbCnf& cnfModifier, \
                                 const int& NConfigs, \
                                 const int& NBars,  \
                                 const string& POT, \
                                 const double& LC, \
                                 const vector<int>& dupFactors, \
                                 vector<string>& elems, \
                                 const vector<int>& nums) {

  /*
   * the idea is NConfigs be iterating over possible combinations for each
   * elements in range
   * and NBars is still how many possible pairs for each configuration
   */
  /* in uniform subMode, input nums is the upper limits of each element */
  int quotient = NConfigs / nProcs;
  int remainder = NConfigs % nProcs;
  int nCycle = remainder ? (quotient + 1) : quotient;
  for (int j = 0; j < nCycle; ++j) {
    for (int i = (j * nProcs); i < ((j + 1) * nProcs); ++i) {
      if ((i % nProcs != me) || (i >= NConfigs)) continue;
      Config c0 = cnfModifier.getFCCConv(LC, elems[0], dupFactors);
      /* get rand ints */
      vector<int> numsVec;
      for (int k = 0; k < nums.size(); ++k) {
        if (elems[k] == "Al") continue;
        int offset;
        if (elems[k] == "X") {
          offset = 1;
        } else {
          offset = 0;
        }
        numsVec.push_back(myRandInt(offset, nums[k]));
      }
      int others = 0; //other than Al
      for (const auto val : numsVec) {
        others += val;
      }
      auto it = numsVec.insert(numsVec.begin(), (c0.natoms - others));

      /* get rand ints end */
      std::default_random_engine rng(std::random_device{}() + me);

      cnfModifier.getRandConfUniformDist(c0, rng, elems, numsVec);
      Config c0copy = c0; //because write POS will sort c0, hence change index

      string baseDir = "config" + to_string(i) + "/s";
      string mkBaseDir = "mkdir -p " + baseDir;
      const char *cmkBaseDir = mkBaseDir.c_str();
      const int dir_err = std::system(cmkBaseDir);
      if (-1 == dir_err) {
        cout << "Error creating directory!\n";
        exit(1);
      }

      cnfModifier.writeCfgData(c0, "config" + to_string(i) + "/s/start.cfg");

      // add small perturbation to break perfect fcc symmetry
      // this method is about to increase the chance
      // to find lower ground states for VASP software
      cnfModifier.cnvprl2pst(c0);
      cnfModifier.perturb(c0);

      map<string, int> elemName = cnfModifier.writePOSCAR(c0, \
                                      "config" + to_string(i) + "/s/POSCAR");
      prepVASPFiles(baseDir, dupFactors, elemName, POT);
      vector<pair<int, int>> pairs = cnfModifier.getPairToSwap(c0copy);

      std::shuffle(std::begin(pairs), std::end(pairs), rng);

      int end = MIN(NBars, pairs.size());
      for (unsigned int k = 0; k < end; ++k) {

        string subDir = "config" + to_string(i) + "/e_" + to_string(k);
        const char *csubDir = subDir.c_str();
        const int dir_err = mkdir(csubDir, \
                                  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (-1 == dir_err) {
          cout << "Error creating directory!\n";
          exit(1);
        }

        Config c1 = cnfModifier.swapPair(c0copy, pairs[k]);
        string name1 = "config" + to_string(i) + "/e_" + to_string(k) + "/";
        cnfModifier.writeCfgData(c1, name1 + "end.cfg");

        // add small perturbation to break perfect fcc symmetry
        // this method is about to increase the chance
        // to find lower ground states for VASP software
        cnfModifier.cnvprl2pst(c1);
        cnfModifier.perturb(c1);

        map<string, int> elemName = cnfModifier.writePOSCAR(c1, \
                                                          name1 + "POSCAR");

        ofstream ofs("log.txt", std::ofstream::app);

        ofs << "config " << i << " end " << k \
            << " pair: " << pairs[k].first \
            << " " << pairs[k].second << "\n";

        ofs.close();

        prepVASPFiles(name1, dupFactors, elemName, POT);
      }
    }
  }
}

void KNHome::createRandomSpecific(gbCnf& cnfModifier, \
                                  const int& NConfigs, \
                                  const int& NBars,  \
                                  const string& POT, \
                                  const double& LC, \
                                  const vector<int>& dupFactors, \
                                  vector<string>& elems, \
                                  const vector<int>& nums) {

  int quotient = NConfigs / nProcs;
  int remainder = NConfigs % nProcs;
  int nCycle = remainder ? (quotient + 1) : quotient;
  for (int j = 0; j < nCycle; ++j) {
    for (int i = (j * nProcs); i < ((j + 1) * nProcs); ++i) {
      if ((i % nProcs != me) || (i >= NConfigs)) continue;
      Config c0 = cnfModifier.getFCCConv(LC, elems[0], dupFactors);

      std::default_random_engine rng(std::random_device{}() + me);
      cnfModifier.getRandConfUniformDist(c0, rng, elems, nums);
      Config c0copy = c0; //because write POS will sort c0, hence change index

      string baseDir = "config" + to_string(i) + "/s";
      string mkBaseDir = "mkdir -p " + baseDir;
      const char *cmkBaseDir = mkBaseDir.c_str();
      const int dir_err = std::system(cmkBaseDir);
      if (-1 == dir_err) {
        cout << "Error creating directory!\n";
        exit(1);
      }

      cnfModifier.writeCfgData(c0, "config" + to_string(i) + \
                                   "/s/start.cfg");

      // add small perturbation to break perfect fcc symmetry
      // this method is about to increase the chance
      // to find lower ground states for VASP software
      cnfModifier.cnvprl2pst(c0);
      cnfModifier.perturb(c0);

      map<string, int> elemName = cnfModifier.writePOSCAR(c0, \
                                      "config" + to_string(i) + "/s/POSCAR");
      prepVASPFiles(baseDir, dupFactors, elemName, POT);
      vector<pair<int, int>> pairs = cnfModifier.getPairToSwap(c0copy);

      std::shuffle(std::begin(pairs), std::end(pairs), rng);

#ifdef DEBUG
      if (me == 0) {
        cout << "config " << i << "\n";
        for (const auto& p : pairs) {
          cout << p.first << " " << p.second << "\n";
        }
      }
#endif

      int end = MIN(NBars, pairs.size());
      for (unsigned int k = 0; k < end; ++k) {

        string subDir = "config" + to_string(i) + "/e_" + to_string(k);
        const char *csubDir = subDir.c_str();
        const int dir_err = mkdir(csubDir, \
                                  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (-1 == dir_err) {
          cout << "Error creating directory!\n";
          exit(1);
        }

        Config c1 = cnfModifier.swapPair(c0copy, pairs[k]);
        string name1 = "config" + to_string(i) + "/e_" + to_string(k) + "/";
        cnfModifier.writeCfgData(c1, name1 + "end.cfg");

        // add small perturbation to break perfect fcc symmetry
        // this method is about to increase the chance
        // to find lower ground states for VASP software
        cnfModifier.cnvprl2pst(c1);
        cnfModifier.perturb(c1);

        map<string, int> elemName = cnfModifier.writePOSCAR(c1, \
                                                          name1 + "POSCAR");

        ofstream ofs("log.txt", std::ofstream::app);

        ofs << "config " << i << " end " << k \
            << " pair: " << pairs[k].first \
            << " " << pairs[k].second << "\n";

        ofs.close();

        prepVASPFiles(name1, dupFactors, elemName, POT);
      }
    }
  }

}
