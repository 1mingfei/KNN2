#include "gbCnf.h"
#include "KNHome.h"
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

vector<pair<int, int>> gbCnf::getPairToSwap(Config& cnf) {
  vector<pair<int, int>> res;
  getNBL(cnf, 3.0);
  for (unsigned int i = 0; i < cnf.vacList.size(); ++i) {
    for (unsigned int j = 0; j < cnf.atoms[cnf.vacList[i]].NBL.size(); ++j) {
      if (cnf.atoms[cnf.vacList[i]].NBL[j] == -1)
        continue;
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
  /* initialize atom type array */
  vector<int> TPArr(cnf.natoms, 0);
  /* start */
  vector<string>::iterator it = std::find(elems.begin(), elems.end(), "X");
  int index = std::distance(elems.begin(), it);
  int nVac = nums[index];
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
    for (int j = 0; j < nums[i]; ++j) {
      if (nblType.size() == NBLSize)
        break;
      nblType.push_back(i);
      --carryOver;
    }
    for (int j = 0; j < carryOver; ++j) {
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

void KNHome::createPreNEB(gbCnf& cnfModifier) {
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
    std::cout << "generating inital and final structures...\n";

  if (subMode == "ordered_cluster") {
    createOrdered(cnfModifier, dupFactors, LC, POT);
  } else if (subMode == "ordered_cluster_random") {
    createOrderedRandom(cnfModifier, dupFactors, LC, \
                        POT, iparams["numDataset"]);
  } else if (subMode == "ordered_cluster_diffcon") {
    createOrderedDiffCon(cnfModifier, dupFactors, LC, \
                        POT, iparams["numDataset"]);
  } else if (subMode == "ordered_cluster_antiphase") {
    createOrderedAntiPhase(cnfModifier, dupFactors, LC, \
                         POT, iparams["numDataset"]);
  } else if (subMode == "random") {
    createRandom(cnfModifier, NConfigs, NBars, POT, LC, \
                 dupFactors, elems, nums);
  } else if (subMode == "uniform") {
    createRandomUniform(cnfModifier, NConfigs, NBars, POT, LC, \
                        dupFactors, elems, nums);
  } else if (subMode == "specific") {
    createRandomSpecific(cnfModifier, NConfigs, NBars, POT, LC, \
                         dupFactors, elems, nums);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (me == 0)
    std::cout << "done generating.\n";
}
