#include "gbCnf.h"
#include "KNHome.h"

#define KB 8.6173303e-5

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

  string fname = sparams["initconfig"];
  RCut = dparams["RCut"];
  maxIter = iparams["maxIter"];
  ntally = iparams["ntally"];
  temperature = dparams["temperature"];
  srand(iparams["randSeed"]);
  prefix = dparams["prefix"];
  c0 = cnfModifier.readCfg(fname);

  /* get initial vacancy position in atomList */
  getVacList();

  /* get neibour list of atoms */
  cnfModifier.wrapAtomPrl(c0);
  cnfModifier.getNBL(c0, 4.5);

  for (auto&& i : vacList) {
    vector<int> tmpVector;
    for (const auto& j : c0.atoms[i].NBL) {
      if (c0.atoms[j].tp == "X")
        continue;
      double dist = cnfModifier.calDistPrl(c0.length, \
                                           c0.atoms[i], \
                                           c0.atoms[j]);
      if (dist <= RCut)
        tmpVector.push_back(j);
    }
    jumpList[i] = tmpVector;

    /* update time */
    time = 0.0;
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
  string modelFname = sparams["kerasModel"];
  k2pModel = Model::load(modelFname);
}

// typedef struct cmp {
//   bool operator() (KMCEvent a, double value) {
//     return (a.getcProb() < value);
//   }
// } comp;

KMCEvent KNHome::selectEvent() {
  double randVal = (double) rand() / (RAND_MAX);
  auto it = mylower_bound(eventList.begin(), eventList.end(), randVal, \
                          [] (KMCEvent a, double value) \
                            {return (a.getcProb() < value);}); 

  if(it == eventList.cend())
    return eventList.back();

// #ifdef DEBUG
  cout << "step " << (step + 1) << " time " << time
       << " prob: " << randVal << " event: " << distance(eventList.begin(), it)\
       << " event cprob: " << it->getcProb() << " jumpPair: " \
       << it->getJumpPair().first << " " << it->getJumpPair().second \
       << endl;
// #endif
  return *it;
}


double KNHome::calRate(Config& c0, \
                       const double& T, \
                       gbCnf& cnfModifier, \
                       pair<int, int> jumpPair) {

  int first = jumpPair.first;
  int second = jumpPair.second;

  if (c0.atoms[first].tp == c0.atoms[second].tp)
    return 0.0;

  vector<string> codes; // atom location in original atom list
  //4.5 for 2NN encoding needed
  vector<vector<string>> encodes = cnfModifier.encodeConfig(c0, \
                      {first, second}, \
                      4.5, \
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

#ifdef DEBUG
    cout << "predicted energies of pair " << first << " " \
         << second << " : ";
#endif

  int nRow = input.size(); // encodings for one jump pair considering symmetry
  int nCol = nRow ? input[0].size() : 0;
  Tensor in{ nRow, nCol };
  for (int i = 0; i < nRow; ++i)
    for (int j = 0; j < nCol; ++j)
      in.data_[i * nCol + j] = input[i][j];

  Tensor out = k2pModel(in);
  double deltaE = 0.0;
  for (int i = 0; i < nRow; ++i) {

#ifdef DEBUG
    cout << std::setprecision(8) << out(i, 0) << " ";
#endif

    deltaE += static_cast<double>(out(i, 0));
  }
  deltaE /= static_cast<double>(nRow);

#ifdef DEBUG
  cout << std::setprecision(8) << deltaE << endl;
#endif

  // double deltaE = (double) rand() / (RAND_MAX);
  return exp(-deltaE / KB / T);
}

void KNHome::buildEventList(gbCnf& cnfModifier) {
  /* build event list */
  // cnfModifier.getNBL(c0, 4.5);

  eventList.clear();
  double sum = 0.0;
  for (int i = 0; i < vacList.size(); ++i) {
    for (int j = 0; j < jumpList[vacList[i]].size(); ++j) {
      /* skip Vac jump to Vac in event list */
      if (c0.atoms[vacList[i]].tp == c0.atoms[jumpList[vacList[i]][j]].tp)
        continue;

      KMCEvent event(make_pair(vacList[i], jumpList[vacList[i]][j]));


      // event.calRate(c0, temperature, RCut, cnfModifier);
      double currRate = calRate(c0, \
                                temperature, \
                                cnfModifier, \
                                make_pair(vacList[i], jumpList[vacList[i]][j]));
      event.setRate(currRate);
      sum += event.getRate();
      eventList.push_back(event);
    }
  }

  /* calculate relative and cumulative probability */
  double curr = 0.0;
  for (int i = 0; i < eventList.size(); ++i) {
    auto&& event = eventList[i];
    event.calProb(sum);
    curr += event.getProb();
    event.setcProb(curr);
  }

  /* update time lapse */
  double tau = prefix * sum * (-log(rand() / static_cast<double>(RAND_MAX)));
  time += tau;

#ifdef DEBUG
  for (int i = 0; i < eventList.size(); ++i) {
    const auto& event = eventList[i];
    cout << setprecision(16) << i << " rate: " << event.getRate() \
         << " cumulative prob: " << event.getcProb() << endl;
  }
  cout << endl;
#endif
}

void KNHome::KMCSimulation(gbCnf& cnfModifier) {
  KMCInit(cnfModifier);
  step = 0;
  cnfModifier.writeCfgData(c0, to_string(step) + ".cfg");

  while (step < maxIter) {
    buildEventList(cnfModifier);
    auto&& event = selectEvent();
    event.exeEvent(c0, jumpList, RCut);
    ++step;
    if (step % ntally == 0)
      cnfModifier.writeCfgData(c0, to_string(step) + ".cfg");

  }
}