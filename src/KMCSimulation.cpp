#include "KNHome.h"
#include "LRUCache.h"
#include "gbCnf.h"

#include <chrono>
using namespace std::chrono;

#define KB 8.6173303e-5
#define KB_INV 11604.5221105
#define NEI_NUMBER 12
#define KEY_SIZE 27

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
    else {
      count = step;
    }
  }
  return first;
}

void KNHome::buildEmbedding() {
  elems = vsparams["elems"];
  switchUnknown = false;
  int base = 0;
  int j = 1;
  for (int i = 1; i <= elems.size(); ++i) {
    if (elems[i - 1] == "Al")
      base = 1;
    if (elems[i - 1] == "Xe") { // unknown element
      embedding[elems[i - 1]] = static_cast<double>(base);
      switchUnknown = true;
    }
    else {
      embedding[elems[i - 1]] = static_cast<double>(j++);
    }
  }

#ifdef DEBUG
  for (int i = 0; i < elems.size(); ++i) {
    cout << embedding[elems[i]] << " ";
  }
  cout << endl;
#endif

  if (switchUnknown) {
    elemsEffectOffset = vdparams["elemsEffectOffset"];
  }

#ifdef DEBUG
  for (int i = 0; i < elems.size(); ++i) {
    cout << elemsEffectOffset[i] << " ";
  }
  cout << endl;
#endif
}

void KNHome::getVacList() {
  for (int i = 0; i < c0.atoms.size(); ++i) {
    if (c0.atoms[i].tp == "X")
      vacList.push_back(i);
  }
}

void KNHome::KMCInit(gbCnf& cnfModifier) {

  buildEmbedding();
  if (me == 0)
    ofs.open("log.txt", std::ofstream::out | std::ofstream::app);

  E_tot = 0.0;
  string fname = sparams["initconfig"];
  RCut = (dparams["RCut"] == 0.0) ? 3.0 : dparams["RCut"];
  RCut2 = 1.65 * RCut;
  maxIter = iparams["maxIter"];
  nTallyConf = (iparams["nTallyConf"] == 0) ? 1000 : iparams["nTallyConf"];
  nTallyOutput = \
           (iparams["nTallyOutput"] == 0) ? 1000 : iparams["nTallyOutput"];
  EDiff = sparams["EDiff"];
  ECutoff = dparams["ECutoff"];

  trapStep = 0;
  temperature = dparams["temperature"];
  srand(iparams["randSeed"]);
  prefix = dparams["prefix"];
  switchLSKMC = false;
  switchLSKMC = bparams["switchLSKMC"];

  c0 = cnfModifier.readCfg(fname);

  /* get initial vacancy position in atomList */
  getVacList();

  /* get neibour list of atoms */
  cnfModifier.wrapAtomPrl(c0);

  if (nProcs == 1) {
    cnfModifier.getNBL_serial(c0, RCut2);
  } else {
    cnfModifier.getNBL(c0, RCut2);
  }

#ifdef DEBUG
  if (me == 0)
    for (int i = 0; i < vacList.size(); ++i) {
      cout << vacList[i] << " size: " << c0.atoms[vacList[i]].FNNL.size() \
      << endl;
      for (int j = 0; j < c0.atoms[vacList[i]].FNNL.size(); ++j)
        cout << c0.atoms[vacList[i]].FNNL[j] << " ";
      cout << endl;
    }
#endif

  /* reading the model binary file, initialize the model */
  string modelFname = sparams["kerasModelBarrier"];
  k2pModelB = Model::load(modelFname);
  if (EDiff == "model") {
    modelFname = sparams["kerasModelEDiff"];
    k2pModelD = Model::load(modelFname);
  }

  /* initialize parameters */
  iter = 0;
  if (sparams["method"] == "restart") {

    time = dparams["startingTime"];
    step = iparams["startingStep"];
    E_tot = dparams["startingEnergy"];
    if (me == 0)
      ofs << "#restarting from step " << step << "\n";

  } else {
    time = (dparams["startingTime"] == 0.0) ? 0.0 : dparams["startingTime"];
    E_tot = \
          (dparams["startingEnergy"] == 0.0) ? 0.0 : dparams["startingEnergy"];
    step = (iparams["startingStep"] == 0) ? 0 : iparams["startingStep"];
    if (me == 0)
      cnfModifier.writeCfgData(c0, to_string(step) + ".cfg");
  }
  if (me == 0)
    ofs << "#step     time     Ediff\n";
    // ofs << "#step     time     Ediff     cachTimes\n";

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

  if (it == eventList.cend()) {
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

double KNHome::updateTime() {
  /* update time elapsed */
  double tau = log(rand() / static_cast<double>(RAND_MAX)) \
               / (prefix * kTot);

  time -= tau;
  return -tau;
}

void KNHome::updateEnergy(const int& eventID) {
  E_tot += eventList[eventID].getEnergyChange();
}

void KNHome::buildEventList_serial(gbCnf& cnfModifier) {
  eventList.clear();
  kTot = 0.0;
  int i, j;

  for (i = 0; i < vacList.size(); ++i) {
    for (j = 0; j < c0.atoms[vacList[i]].FNNL.size(); ++j) {

      int iFirst = vacList[i];
      int iSecond = c0.atoms[vacList[i]].FNNL[j];

      KMCEvent event(make_pair(iFirst, iSecond));

      vector<double> currBarrier;

      if (lru->getSize()) {
        currBarrier = cnfModifier.calBarrierAndEdiff_LRU(c0, \
                                  temperature, \
                                  RCut2, \
                                  EDiff, \
                                  embedding, \
                                  k2pModelB, \
                                  k2pModelD, \
                                  make_pair(iFirst, iSecond), \
                                  switchUnknown, \
                                  elems, \
                                  elemsEffectOffset, \
                                  lru);
      } else {
        currBarrier = cnfModifier.calBarrierAndEdiff(c0, \
                                  temperature, \
                                  RCut2, \
                                  EDiff, \
                                  embedding, \
                                  k2pModelB, \
                                  k2pModelD, \
                                  make_pair(iFirst, iSecond), \
                                  switchUnknown, \
                                  elems, \
                                  elemsEffectOffset);
      }

      if (c0.atoms[iFirst].tp == c0.atoms[iSecond].tp) {
        event.setRate(0.0);
        event.setEnergyChange(0.0);
        event.setBarrier(0.0);
      } else {
        event.setRate(exp(-currBarrier[0] * KB_INV / temperature));
        event.setEnergyChange(currBarrier[1]);
        event.setBarrier(currBarrier[0]);
      }

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

void KNHome::buildEventList(gbCnf& cnfModifier) {

  eventList.clear();
  eventList.resize(NEI_NUMBER);
  kTot = 0.0;

  int iFirst = vacList[0];

  int NI = c0.atoms[iFirst].FNNL.size();

  int quotient = NI / nProcs;
  int remainder = NI % nProcs;
  int nCycle = remainder ? (quotient + 1) : quotient;

  double** data;
  data = alloc_2d_array<double> (nCycle * nProcs, 5);

  // smallest buff for gathering in each cycle
  double* buffData = new double [5];

#ifdef DEBUG_MPI
  if (me == 0) {
    for (int i = 0; i < NI; ++i) {
      cout << c0.atoms[iFirst].FNNL[i] << " " \
           << c0.atoms[c0.atoms[iFirst].FNNL[i]].tp << " ";
    }
    cout << "\n";
  }
#endif


  for (int j = 0; j < nCycle; ++j) {

    for (int i = (j * nProcs); i < ((j + 1) * nProcs); ++i) {

      int iSecond = c0.atoms[iFirst].FNNL[i];

      if ((me == 0) && (i % nProcs != 0)) {
        MPI_Recv(&data[i][0], 5, MPI_DOUBLE, (i % nProcs), 0, \
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }

      if (i % nProcs != me) continue;

      for (int k = 0; k < 5; ++k)
        buffData[j] = 0.0;

      if (i >= NI) continue;

#ifdef DEBUG_MPI
      cout << "processor #" << me << " iteration: "
           << iter << " neighbor #" << i << " pair: "
           << iFirst << " " << iSecond << " " << c0.atoms[iSecond].tp << "\n";
#endif

      vector<double> currBarrier;

      if (lru->getSize()) {
        currBarrier = cnfModifier.calBarrierAndEdiff_LRU(c0, \
                                  temperature, \
                                  RCut2, \
                                  EDiff, \
                                  embedding, \
                                  k2pModelB, \
                                  k2pModelD, \
                                  make_pair(iFirst, iSecond), \
                                  switchUnknown, \
                                  elems, \
                                  elemsEffectOffset, \
                                  lru);
      } else {
        currBarrier = cnfModifier.calBarrierAndEdiff(c0, \
                                  temperature, \
                                  RCut2, \
                                  EDiff, \
                                  embedding, \
                                  k2pModelB, \
                                  k2pModelD, \
                                  make_pair(iFirst, iSecond), \
                                  switchUnknown, \
                                  elems, \
                                  elemsEffectOffset);
      }

#ifdef DEBUG_MPI
      cout << "processor #" << me << " iteration: " \
           << iter << " neighbor #" << i \
           << " calulate barrier done.\n";
#endif

      if (c0.atoms[iFirst].tp != c0.atoms[iSecond].tp) {
        buffData[0] = exp(-currBarrier[0] * KB_INV / temperature); // rate
        buffData[1] = currBarrier[1];                              // Ediff
        buffData[2] = currBarrier[0];                              // Ebarr
        buffData[3] = static_cast<double>(iFirst);                 // jmpPair1
        buffData[4] = static_cast<double>(iSecond);                // jmpPair2
      }

      if (me != 0) {
        MPI_Send(&buffData[0], 5, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      } else {
        for (int ii = 0; ii < 5; ++ii)
          data[j * nProcs + i % nProcs][ii] = buffData[ii];
      }

    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  MPI_Bcast(&data[0][0], nCycle * nProcs * 5, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  for (int i = 0; i < NEI_NUMBER; ++i) {
    kTot += data[i][0];
    eventList[i].setRate(data[i][0]);
    eventList[i].setEnergyChange(data[i][1]);
    eventList[i].setBarrier(data[i][2]);
    eventList[i].setJumpPair(static_cast<int>(data[i][3]), \
                             static_cast<int>(data[i][4]));
  }

  /* calculate relative and cumulative probability */
  double curr = 0.0;
  for (int i = 0; i < eventList.size(); ++i) {
    auto&& event = eventList[i];
    event.calProb(kTot);
    curr += event.getProb();
    event.setcProb(curr);
  }

  // MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG_MPI
  cout << "processor #" << me << " iteration: "
       << iter << " set eventList done.\n";
#endif

  /*free smallest buffer*/
  delete [] buffData;

  /*free largest data set*/
  free(data[0]);
  free(data);

}

void KNHome::KMCSimulation(gbCnf& cnfModifier) {

  KMCInit(cnfModifier);

  while (iter < maxIter) {

    if (nProcs == 1) {
      buildEventList_serial(cnfModifier);
    } else {
      buildEventList(cnfModifier);
    }

#ifdef DEBUG_MPI
    cout << "built event list processor #" << me << " iteration: "
         << iter << " eventList size: " << eventList.size() << "\n";
#endif

    int eventID = 0;

    if (me == 0) {
      auto&& event = selectEvent(eventID);
    }

    MPI_Bcast(&eventID, 1, MPI_INT, 0, MPI_COMM_WORLD);

#ifdef DEBUG_MPI
    cout << "iteration: " << iter << "processor #" << me \
         << " selected event #" << eventID << "\n";
#endif

    eventList[eventID].exeEvent(c0, RCut); // event updated

#ifdef DEBUG_MPI
    cout << "iteration: " << iter << "processor #" << me \
         << " event #" << eventID << " jump pair " \
         << eventList[eventID].getJumpPair().first << "-" \
         << eventList[eventID].getJumpPair().second << "\n";
#endif

    if (me == 0) {
      double oneStepTime = updateTime();
      updateEnergy(eventID);
    }

    ++step;
    ++iter;

    if (me == 0) {
      if (step % nTallyOutput == 0)
        ofs << std::setprecision(7) << step << " " << time << " " \
            << E_tot << " " << lru->getCt() <<endl;

      if (step % nTallyConf == 0)
        cnfModifier.writeCfgData(c0, to_string(step) + ".cfg");
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

}


vector<double> gbCnf::calBarrierAndEdiff(Config& c0, \
                              const double& T, \
                              const double& RCut2, \
                              const string& EDiff, \
                              unordered_map<string, double>& embedding, \
                              Model& k2pModelB, \
                              Model& k2pModelD, \
                              const pair<int, int>& jumpPair, \
                              const bool& switchUnknown, \
                              vector<string>& elems, \
                              const vector<double>& elemsEffectOffset) {

  int first = jumpPair.first;
  int second = jumpPair.second;

  vector<string> codes; // atom location in original atom list
  vector<vector<string>> encodes = encodeConfig(c0, \
                                                {first, second}, \
                                                RCut2, \
                                                codes, \
                                                {first, second}, \
                                                false);
  vector<vector<double>> input(encodes.size(), \
                               vector<double>(encodes[0].size(), 0.0));

  double Eactivate = 0.0;
  double EactivateBack = 0.0;

  for (int i = 0; i < encodes.size(); ++i)
    for (int j = 0; j < encodes[i].size(); ++j)
      input[i][j] = embedding[encodes[i][j]];

  int nRow = input.size(); // encodings for one jump pair considering symmetry
  Tensor in{ nRow, KEY_SIZE };
  Tensor inBack{ nRow, KEY_SIZE };

  for (int i = 0; i < nRow; ++i) {
    inBack.data_[i * KEY_SIZE] = input[i][0];
    for (int j = 0; j < KEY_SIZE; ++j) {
      in.data_[i * KEY_SIZE + j] = input[i][j];
      if (j == 0)
        continue;
      inBack.data_[i * KEY_SIZE + j] = input[i][KEY_SIZE - j];
    }
  }

  Tensor outB = k2pModelB(in);
  Tensor outBBack = k2pModelB(inBack);

  for (int i = 0; i < nRow; ++i) {
    Eactivate += static_cast<double>(outB(i, 0));
    EactivateBack += static_cast<double>(outBBack(i, 0));
  }

  Eactivate /= static_cast<double>(nRow);
  EactivateBack /= static_cast<double>(nRow);

  double Ediff = Eactivate - EactivateBack;
  if (switchUnknown)  {
    Eactivate += offsetBarrier(c0, elems, elemsEffectOffset, {first, second});
  }
  return {Eactivate, Ediff};
}

vector<double> gbCnf::calBarrierAndEdiff_LRU(Config& c0, \
                              const double& T, \
                              const double& RCut2, \
                              const string& EDiff, \
                              unordered_map<string, double>& embedding, \
                              Model& k2pModelB, \
                              Model& k2pModelD, \
                              const pair<int, int>& jumpPair, \
                              const bool& switchUnknown, \
                              vector<string>& elems, \
                              const vector<double>& elemsEffectOffset, \
                              LRUCache* lru) {

  int first = jumpPair.first;
  int second = jumpPair.second;

  vector<string> codes; // atom location in original atom list
  vector<vector<string>> encodes = encodeConfig(c0, \
                                                {first, second}, \
                                                RCut2, \
                                                codes, \
                                                {first, second}, \
                                                false);


  vector<array<int, KEY_SIZE>> input;
  vector<array<int, KEY_SIZE>> inputBack;

  double Eactivate = 0.0;
  double EactivateBack = 0.0;

  int nRow = encodes.size();

  for (int i = 0; i < nRow; ++i) {

    array<int, KEY_SIZE> tmpVec;
    array<int, KEY_SIZE> tmpVecBack;
    tmpVecBack[0] = embedding[encodes[i][0]];

    for (int j = 0; j < KEY_SIZE; ++j) {
      tmpVec[j] = embedding[encodes[i][j]];
      if (j == 0) continue;
      tmpVecBack[j] = embedding[encodes[i][KEY_SIZE - j]];
    }
    if (lru->check(tmpVec)) {
      Eactivate += lru->getBarrier(tmpVec);
    } else {
      input.push_back(tmpVec);
    }

    if (lru->check(tmpVecBack)) {
      EactivateBack += lru->getBarrier(tmpVecBack);
    } else {
      inputBack.push_back(tmpVecBack);
    }
  }

  nRow = input.size(); // encodings for one jump pair considering symmetry
  int nRowBack = inputBack.size();

  if (nRow) {
    Tensor in{nRow, KEY_SIZE};
    for (int i = 0; i < nRow; ++i) {
      for (int j = 0; j < KEY_SIZE; ++j) {
        in.data_[i * KEY_SIZE + j] = input[i][j];
      }
    }

    Tensor outB = k2pModelB(in);

    for (int i = 0; i < nRow; ++i) {
      double tmpEa = static_cast<double>(outB(i, 0));
      Eactivate += tmpEa;
      lru->add(make_pair(input[i], tmpEa));
    }
  }
  if (nRowBack) {
    Tensor inBack{nRowBack, KEY_SIZE};
    for (int i = 0; i < nRowBack; ++i) {
      for (int j = 0; j < KEY_SIZE; ++j) {
        inBack.data_[i * KEY_SIZE + j] = inputBack[i][j];
      }
    }
    Tensor outBBack = k2pModelB(inBack);
    for (int i = 0; i < nRowBack; ++i) {
      double tmpEaBack = static_cast<double>(outBBack(i, 0));
      EactivateBack += tmpEaBack;
      lru->add(make_pair(inputBack[i], tmpEaBack));
    }
  }

  Eactivate /= static_cast<double>(encodes.size());
  EactivateBack /= static_cast<double>(encodes.size());

  double Ediff = Eactivate - EactivateBack;
  if (switchUnknown)  {
    Eactivate += offsetBarrier(c0, elems, elemsEffectOffset, {first, second});
  }
  return {Eactivate, Ediff};
}

double gbCnf::offsetBarrier(const Config& c0, \
                            vector<string>& elems, \
                            const vector<double>& elemsEffectOffset, \
                            const pair<int, int>& jumpPair) {

  const int iFirst = jumpPair.first;
  const int iSecond = jumpPair.second;

  double res = 0.0;

  if (c0.atoms[iSecond].tp == "Xe") {
    // This map should be the bonding numbers of X-Al, X-Mg, X-Zn,
    // key is Al, Mg or Zn
    unordered_map<string, int> mp;
    for (int i = 0; i < NEI_NUMBER; ++i) {
      int n = c0.atoms[iFirst].FNNL[i];
      --mp[c0.atoms[n].tp];

      int m = c0.atoms[iSecond].FNNL[i];
      ++mp[c0.atoms[m].tp];
    }

#ifdef DEBUG_OFFSET
    cout << "jump element is pseudo\n";
#endif

    for (int i = 0; i < elems.size(); ++i) {
#ifdef DEBUG_OFFSET
      cout << "#debug here #" << mp[elems[i]] \
           <<  " " << elemsEffectOffset[i] << endl;
#endif
      res += (static_cast<double>(mp[elems[i]]) * elemsEffectOffset[i]);
    }
  } else {
    // This map should be the bonding numbers of Al-X or Mg-X or Zn-X
    // then only need to count how many ?-X bond change before and after
    int count = 0;
    for (int i = 0; i < NEI_NUMBER; ++i) {
      int n = c0.atoms[iFirst].FNNL[i];
      if (c0.atoms[n].tp == "Xe")
        --count;

      int m = c0.atoms[iSecond].FNNL[i];
      if (c0.atoms[m].tp == "Xe")
        ++count;
    }
    vector<string>::iterator it = std::find(elems.begin(), elems.end(), \
                                            c0.atoms[iSecond].tp);
    int index = std::distance(elems.begin(), it);

    res += (static_cast<double>(count) * elemsEffectOffset[index]);
#ifdef DEBUG_OFFSET
    cout << "jump element is pseudo\n";
    cout << "#debug here #" << count \
         <<  " " << elemsEffectOffset[index] << endl;
#endif
  }

#ifdef DEBUG_OFFSET
  cout << "#debug here #" << res << endl;
#endif

  return res;
}
