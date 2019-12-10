#include "gbCnf.h"
#include "KNHome.h"

#define KB 8.6173303e-5
#define NEI_NUMBER 12

/*  first element in the range [first, last)
 *  that is not less than (i.e. greater or equal to) value */
template<class ForwardIt, class T, class Compare>
inline ForwardIt mylower_bound(ForwardIt first, \
                               ForwardIt last, \
                               const T& value,\
                               Compare comp) {
  ForwardIt it;
  typename std::iterator_traits<ForwardIt>::difference_type count, step;
  count = std::distance(first, last);

  while (count > 0) {
    it = first;
    step = count / 2;
    std::advance(it, step);
    if (comp(*it, value)) {
      first = ++it;
      count -= step + 1;
    }
    else
      count = step;
  }
  return first;
}

void KNHome::buildEmbedding() {
  vector<string> elems = vsparams["elems"];
  for (int i = 1; i <= elems.size(); ++i)
    embedding[elems[i - 1]] = static_cast<double>(i);
}



void KNHome::getVacList() {
  for (int i = 0; i < c0.atoms.size(); ++i) {
    if (c0.atoms[i].tp == "X")
      vacList.push_back(i);
  }
}

void KNHome::KMCInit(gbCnf& cnfModifier) {

  buildEmbedding();
  E_tot = 0.0;
  string fname = sparams["initconfig"];
  RCut = dparams["RCut"];
  RCut2 = 1.65 * RCut;
  maxIter = iparams["maxIter"];
  nTallyConf = iparams["nTallyConf"];
  nTallyOutput = iparams["nTallyOutput"];
  switchEngy = bparams["EDiff"];

  temperature = dparams["temperature"];
  srand(iparams["randSeed"]);
  prefix = dparams["prefix"];
  c0 = cnfModifier.readCfg(fname);

  /* get initial vacancy position in atomList */
  getVacList();

  /* get neibour list of atoms */
  cnfModifier.wrapAtomPrl(c0);
  cnfModifier.getNBL(c0, RCut2);

  for (auto&& i : vacList) {
    vector<int> tmpVector;
    for (const auto& j : c0.atoms[i].NBL) {
      // if (c0.atoms[j].tp == "X")
      //   continue;
      double dist = cnfModifier.calDistPrl(c0.length, \
                                           c0.atoms[i], \
                                           c0.atoms[j]);
      if (dist <= RCut)
        tmpVector.push_back(j);
    }
    jumpList[i] = tmpVector;
  }

#ifdef DEBUG
  for (int i = 0; i < vacList.size(); ++i) {
    cout << vacList[i] << " size: " << jumpList[vacList[i]].size() << endl;
    for (int j = 0; j < jumpList[vacList[i]].size(); ++j)
      cout << jumpList[vacList[i]][j] << " ";
    cout << endl;
  }
#endif

  /* reading the model binary file, initialize the model */
  string modelFname = sparams["kerasModelBarrier"];
  k2pModelB = Model::load(modelFname);
  if (switchEngy) {
    modelFname = sparams["kerasModelEDiff"];
    k2pModelD = Model::load(modelFname);
  }
  /* initialize time */
  if (sparams["method"] == "restart") {

    time = dparams["startingTime"];
    step = iparams["startingStep"];
    iter = 0;
    cout << "#restarting from step " << step << "\n";

  } else {
    iter = 0;
    time = 0.0;
    step = 0;
  }

  cnfModifier.writeCfgData(c0, to_string(step) + ".cfg");
  cout << "#step     time     Ediff     jumpFrom     jumpTo\n";

}

// typedef struct cmp {
//   bool operator() (KMCEvent a, double value) {
//     return (a.getcProb() < value);
//   }
// } comp;

KMCEvent KNHome::selectEvent(int& dist) {
  double randVal = (double) rand() / (RAND_MAX);
  auto it = mylower_bound(eventList.begin(), eventList.end(), randVal, \
                          [] (KMCEvent a, double value) \
                            {return (a.getcProb() < value);});

  if(it == eventList.cend()) {
    dist = eventList.size() - 1;
    return eventList.back();
  }

  dist = distance(eventList.begin(), it);
#ifdef DEBUGJUMP
  cout << "step " << (step + 1) << " time " << time
       << " prob: " << randVal << " event: " << dist \
       << " event cprob: " << it->getcProb() << " jumpPair: " \
       << it->getJumpPair().first << " " << it->getJumpPair().second \
       << endl;
#endif
  return *it;
}


vector<double> KNHome::calRate(Config& c0, \
                               const double& T, \
                               gbCnf& cnfModifier, \
                               pair<int, int> jumpPair) {

  int first = jumpPair.first;
  int second = jumpPair.second;

  if (c0.atoms[first].tp == c0.atoms[second].tp)
    return {0.0, 0.0};

  vector<string> codes; // atom location in original atom list
  //RCut2 for 2NN encoding needed
  vector<vector<string>> encodes = cnfModifier.encodeConfig(c0, \
                                                            {first, second}, \
                                                            RCut2, \
                                                            codes, \
                                                            {first, second}, \
                                                            false);
  vector<vector<double>> input(encodes.size(), \
                               vector<double>(encodes[0].size(), 0.0));

  for (int i = 0; i < encodes.size(); ++i)
    for (int j = 0; j < encodes[i].size(); ++j)
      input[i][j] = embedding[encodes[i][j]];

#ifdef DEBUG
  for (int i = 0; i < encodes.size(); ++i) {
    for (int j = 0; j < encodes[i].size(); ++j)
      cout << encodes[i][j] << " ";
    cout << endl;
  }
  for (int i = 0; i < input.size(); ++i) {
    for (int j = 0; j < input[i].size(); ++j)
      cout << input[i][j] << " ";
    cout << endl;
  }
#endif

  int nRow = input.size(); // encodings for one jump pair considering symmetry
  int nCol = nRow ? input[0].size() : 0;
  Tensor in{ nRow, nCol };
  for (int i = 0; i < nRow; ++i)
    for (int j = 0; j < nCol; ++j)
      in.data_[i * nCol + j] = input[i][j];
  if (switchEngy) {
    Tensor outB = k2pModelB(in);
    Tensor outD = k2pModelD(in);

    double deltaE = 0.0;
    double tmpEdiff = 0.0;
    for (int i = 0; i < nRow; ++i) {
      deltaE += static_cast<double>(outB(i, 0));
      tmpEdiff += static_cast<double>(outD(i, 0));
    }
    deltaE /= static_cast<double>(nRow);
    tmpEdiff /= static_cast<double>(nRow);

    return {exp(-deltaE / KB / T), tmpEdiff};
  } else {
    Tensor outB = k2pModelB(in);
    // Tensor outD = k2pModelD(in);

    double deltaE = 0.0;
    // double tmpEdiff = 0.0;
    for (int i = 0; i < nRow; ++i) {
      deltaE += static_cast<double>(outB(i, 0));
      // tmpEdiff += static_cast<double>(outD(i, 0));
    }
    deltaE /= static_cast<double>(nRow);
    // tmpEdiff /= static_cast<double>(nRow);

    return {exp(-deltaE / KB / T), 0.0};
  }
}

void KNHome::updateTime() {
  /* update time elapsed */
  double tau = log(rand() / static_cast<double>(RAND_MAX)) \
               / (prefix * kTot);
  time -= tau;
}

void KNHome::updateEnergy(const int& eventID) {
  E_tot += eventList[eventID].getEnergyChange();
}

void KNHome::buildEventList(gbCnf& cnfModifier) {
  eventList.clear();
  kTot = 0.0;
  for (int i = 0; i < vacList.size(); ++i) {
    for (int j = 0; j < jumpList[vacList[i]].size(); ++j) {
      /* skip Vac jump to Vac in event list */

      int iFirst = vacList[i];
      int iSecond = jumpList[vacList[i]][j];

      // if (c0.atoms[iFirst].tp == c0.atoms[iSecond].tp)
      //   continue;

      KMCEvent event(make_pair(iFirst, iSecond));
      string tmpHash = to_string(iFirst) + "_" + to_string(iSecond);
      eventListMap[tmpHash] = i * jumpList[0].size() + j;

      vector<double> currRate = calRate(c0, \
                                temperature, \
                                cnfModifier, \
                                make_pair(iFirst, iSecond));

      event.setRate(currRate[0]);
      event.setEnergyChange(currRate[1]);

      kTot += event.getRate();
      eventList.push_back(event);
    }
  }

  /* calculate relative and cumulative probability */
  double curr = 0.0;
  for (int i = 0; i < eventList.size(); ++i) {
    auto&& event = eventList[i];
    event.calProb(kTot);
    curr += event.getProb();
    event.setcProb(curr);
  }

#ifdef DEBUGJUMP
  for (int i = 0; i < eventList.size(); ++i) {
    const auto& event = eventList[i];
    cout << setprecision(12) << i << " rate: " << event.getRate() \
         << " cumulative prob: " << event.getcProb() << endl;
  }
  cout << endl;
#endif
}

void KNHome::updateEventList(gbCnf& cnfModifier, \
                             const pair<int, int>& jumpPair, \
                             const int& eventID) {
  eventList.clear();
  kTot = 0.0;
  for (int i = 0; i < vacList.size(); ++i) {
    for (int j = 0; j < jumpList[vacList[i]].size(); ++j) {
      /* skip Vac jump to Vac in event list */

      int iFirst = vacList[i];
      int iSecond = jumpList[vacList[i]][j];

      // if (c0.atoms[iFirst].tp == c0.atoms[iSecond].tp)
      //   continue;

      KMCEvent event(make_pair(iFirst, iSecond));
      string tmpHash = to_string(iFirst) + "_" + to_string(iSecond);
      eventListMap[tmpHash] = i * jumpList[0].size() + j;

      vector<double> currRate = calRate(c0, \
                                temperature, \
                                cnfModifier, \
                                make_pair(iFirst, iSecond));

      event.setRate(currRate[0]);
      event.setEnergyChange(currRate[1]);

      kTot += event.getRate();
      eventList.push_back(event);
    }
  }

  /* calculate relative and cumulative probability */
  double curr = 0.0;
  for (int i = 0; i < eventList.size(); ++i) {
    auto&& event = eventList[i];
    event.calProb(kTot);
    curr += event.getProb();
    event.setcProb(curr);
  }

#ifdef DEBUGJUMP
  for (int i = 0; i < eventList.size(); ++i) {
    const auto& event = eventList[i];
    cout << setprecision(12) << i << " rate: " << event.getRate() \
         << " cumulative prob: " << event.getcProb() << endl;
  }
  cout << endl;
#endif

  /* new implementation here */
  /*
      1. remove old keys from hash
      2. calculate where to start and stop, then do it.
      3. update any previous that has iFirst or iSecond that need to be updated.
   */

  // step 1
  // int iFirst = jumpPair.first;
  // int iSecond = jumpPair.second;
  // for (const auto& i : oldJumpList[iFirst]) {
  //   string tmpName = to_string(iFirst) + "_" + to_string(i);
  //   eventListMap.erase(tmpName);
  // }

  //
  // int start = eventID / NEI_NUMBER;
  // for (int i = 0; i < JumpList[iFirst].size(); ++i) {
  //   kTot -= eventList[start + i];
  //   KMCEvent event(make_pair(iFirst, JumpList[iFirst][i]));

  //   string tmpHash = to_string(iFirst) + "_" + to_string(JumpList[iFirst][i]);
  //   eventListMap[tmpHash] = start + i;

  //   double currRate = calRate(c0, \
  //                             temperature, \
  //                             cnfModifier, \
  //                             make_pair(iFirst, JumpList[iFirst][i]));
  //   event.setRate(currRate);
  //   kTot += event.getRate();
  //   eventList[start + i] = event;
  // }

  // for (pair<int, vector<int>>&& elem : jumpList) {
  //   if (elem.first == first)
  //     continue;
  //   for (int j = 0; j < elem.second.size(); ++j) {
  //   }
  // }

}

void KNHome::KMCSimulation(gbCnf& cnfModifier) {
  KMCInit(cnfModifier);

  // buildEventList(cnfModifier);
  while (iter < maxIter) {
    buildEventList(cnfModifier);
    int eventID = 0;
    auto&& event = selectEvent(eventID);
    oldJumpList = jumpList;
    event.exeEvent(c0, jumpList, RCut); // event updated
    updateTime();
    updateEnergy(eventID);
    // updateEventList(cnfModifier, event.getJumpPair(), eventID);
    ++step;
    ++iter;

    if (step % nTallyOutput == 0)
    cout << std::setprecision(7) << step << " " << time << " " \
         << E_tot << " " \
         << event.getJumpPair().first << " " << event.getJumpPair().second \
         << endl;

    if (step % nTallyConf == 0)
      cnfModifier.writeCfgData(c0, to_string(step) + ".cfg");

  }
}