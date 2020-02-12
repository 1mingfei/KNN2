#include "gbCnf.h"
#include "KNHome.h"
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

// random generator function:
inline int myRandInt(int minVal, int maxVal) {
  return minVal + (std::rand() % static_cast<int>(maxVal - minVal + 1));
}

vector<pair<int, int>> gbCnf::getPairToSwap(Config& cnf) {
  vector<pair<int, int>> res;
  getNBL(cnf, 3.0);
  for (unsigned int i = 0; i < cnf.vacList.size(); ++i) {
    for (unsigned int j = 0; j < cnf.atoms[cnf.vacList[i]].NBL.size(); ++j) {
      if (cnf.atoms[j].tp != cnf.atoms[cnf.vacList[i]].tp) {
        res.push_back( std::make_pair(cnf.vacList[i], \
              cnf.atoms[cnf.vacList[i]].NBL[j]) );
      }
    }
  }
  return res;
}

Config gbCnf::swapPair(const Config& c0, pair<int, int> atomPair) {
  int idx0 = atomPair.first;
  int idx1 = atomPair.second;
  Config c1 = c0;

  swap(c1.atoms[idx0].prl, c1.atoms[idx1].prl);
  swap(c1.atoms[idx0].pst, c1.atoms[idx1].pst);
  return c1;
}

void gbCnf::getRandConf(Config& cnf,\
                        std::default_random_engine& rng, \
                        const vector<string>& elems,\
                        const vector<int>& nums) {
  // assert(cnf.natoms == std::accumulate(nums.begin(), nums.end(), 0));
  vector<int> TPArr(nums[0], 0);
  for (unsigned int i = 1; i < nums.size(); ++i) {
    for (unsigned int j = 0; j < nums[i]; ++j) {
      TPArr.push_back(i);
    }
  }

  std::shuffle(std::begin(TPArr), std::end(TPArr), rng);

  for (unsigned int i = 0; i < TPArr.size(); ++i) {
    for (unsigned int j = 0; j < elems.size(); ++j) {
      if (TPArr[i] == j) {
        cnf.atoms[i].tp = elems[j];
      }
    }
  }
  std::sort(cnf.atoms.begin(), cnf.atoms.end());
  for (unsigned int i = 0; i < TPArr.size(); ++i) {
    if (cnf.atoms[i].tp == "X") {
      cnf.vacList.push_back(i);
    }
  }
}

void gbCnf::getRandConfUniformDist(Config& cnf,\
                                   std::default_random_engine& rng, \
                                   vector<string>& elems,\
                                   const vector<int>& nums) {
  // assert(cnf.natoms == std::accumulate(nums.begin(), nums.end(), 0));
  /* initialize atom type array */
  vector<int> TPArr(cnf.natoms, 0);
  //TO-DO: looking for vacancy locations, and corresponding NBL
  //       change NB atom types and countdown till 0
  /* start */
  vector<string>::iterator it = std::find(elems.begin(), elems.end(), "X");
  int index = std::distance(elems.begin(), it);
  int nVac = nums[index];
  // assert(nVac > 0);
  // assign the atom#0 to be the vacancy site
  // and keep assigning neigbhors of it by numbers of Mg Zn and X
  TPArr[0] = index;
  --nVac;

  getNBL(cnf, 5.0);

  int NBLSize = cnf.atoms[0].NBL.size();
  set<int> nblSet;
  for (const auto & i : cnf.atoms[0].NBL) {
    nblSet.insert(i);
  }
  vector<int> nblType;
  vector<int> resType;

  //int locRand = myRandInt(0, nVac);
  int locRand = nVac;
  int carryOver = nVac;

  /* do vacancy separately because -1 vacancy */
  for (unsigned int i = 0; i < locRand; ++i) {
    if (nblType.size() == NBLSize)
      break;
    nblType.push_back(index);
    --carryOver;
  }
  for (int i = 0; i < carryOver; ++i) {
    resType.push_back(index);
  }

  for (unsigned int i = 0; i < elems.size(); ++i) {
    if ((elems[i] == "X") || (elems[i] == "Al"))
      continue;
    carryOver = nums[i];
    //for (int j = 0; j < myRandInt(0, nums[i]); ++j) {
    for (int j = 0; j < nums[i]; ++j) {
      if (nblType.size() == NBLSize)
        break;
      nblType.push_back(i);
      --carryOver;
    }
    for (int j = 0; j < carryOver; ++ j) {
      resType.push_back(i);
    }
  }
  /* make compensation to the lists */
  for (unsigned int i = nblType.size(); i < NBLSize; ++i) {
    nblType.push_back(0);
  }
  for (unsigned int i = resType.size(); i < (cnf.natoms - NBLSize - 1); ++i) {
    resType.push_back(0);
  }
  std::shuffle(std::begin(nblType), std::end(nblType), rng);
  std::shuffle(std::begin(resType), std::end(resType), rng);

  // type alignment (not starting from 0 becasue 0 is the default vacancy)
  int countNBL = 0;
  int countRes = 0;
  for (int i = 1; i < cnf.natoms; ++i) {
    if (nblSet.find(i) != nblSet.end()) {
      TPArr[i] = nblType[countNBL++];
    } else {
      TPArr[i] = resType[countRes++];
    }
  }

  for (unsigned int i = 0; i < TPArr.size(); ++i) {
    for (unsigned int j = 0; j < elems.size(); ++j) {
      if (TPArr[i] == j) {
        cnf.atoms[i].tp = elems[j];
      }
    }
  }

  std::sort(cnf.atoms.begin(), cnf.atoms.end());
  for (unsigned int i = 0; i < TPArr.size(); ++i) {
    if (cnf.atoms[i].tp == "X") {
      cnf.vacList.push_back(i);
    }
  }
}

void KNHome::createPreNEB() {
  gbCnf cnfModifier(sparams);
  vector<int> dupFactors = viparams["factors"];
  double LC = dparams["LC"];
  vector<string> elems = vsparams["elems"];
  vector<int> nums = viparams["nums"];
  int NConfigs = iparams["NConfigs"];
  int NBars = iparams["NBarriers"];
  const string& subMode = sparams["method"];
  const string& POT = sparams["POT"];
  std::default_random_engine rng(std::random_device{}() + me);

  MPI_Barrier(MPI_COMM_WORLD);
  if (me == 0)
    std::cout << "generating inital and final structures\n";
  if (subMode == "ordered_cluster") {
    createOrdered(cnfModifier, dupFactors, LC, POT);
  } else if (subMode == "ordered_cluster_random"){
    createOrderedRandom(cnfModifier, dupFactors, LC, POT);
  } else if (subMode == "random") {
    int quotient = NConfigs / nProcs;
    int remainder = NConfigs % nProcs;
    int nCycle = remainder ? (quotient + 1) : quotient;
    for (int j = 0; j < nCycle; ++j) {
      for (int i = (j * nProcs); i < ((j + 1) * nProcs); ++i) {
        if ((i % nProcs != me) || (i >= NConfigs)) continue;

        Config c0 = cnfModifier.getFCCConv(LC, elems[0], dupFactors);
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

#ifdef DEBUG
        if (me == 0) {
          cout << "config " << i << "\n";
          for (const auto& p : pairs) {
            cout << p.first << " " << p.second << "\n";
          }
        }
#endif

        // auto rng = std::default_random_engine {};
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
  } else if (subMode == "uniform") {
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

#ifdef DEBUG
        if (me == 0) {
          cout << "config " << i << "\n";
          for (const auto& p : pairs) {
            cout << p.first << " " << p.second << "\n";
          }
        }
#endif

        // auto rng = std::default_random_engine {};
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
  } else if (subMode == "specific") {
    int quotient = NConfigs / nProcs;
    int remainder = NConfigs % nProcs;
    int nCycle = remainder ? (quotient + 1) : quotient;
    for (int j = 0; j < nCycle; ++j) {
      for (int i = (j * nProcs); i < ((j + 1) * nProcs); ++i) {
        if ((i % nProcs != me) || (i >= NConfigs)) continue;
        Config c0 = cnfModifier.getFCCConv(LC, elems[0], dupFactors);
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

#ifdef DEBUG
        if (me == 0) {
          cout << "config " << i << "\n";
          for (const auto& p : pairs) {
            cout << p.first << " " << p.second << "\n";
          }
        }
#endif

        // auto rng = std::default_random_engine {};
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
}
