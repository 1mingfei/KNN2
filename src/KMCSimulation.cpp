#include "gbCnf.h"
#include "KNHome.h"

void KNHome::getVacList() {
  for (int i = 0; i < c0.atoms.size(); ++i) {
    if (c0.atoms[i].tp == "X")
      vacList.push_back(i);
  }
}


void KNHome::KMCInit() {

  gbCnf cnfModifier(*this);

  string fname = sparams["initconfig"];
  RCut = dparams["RCut"];
  maxIter = iparams["maxIter"];
  temperature = dparams["temperature"];
  c0 = cnfModifier.readCfg(fname);

  /* get initial vacancy position in atomList */
  getVacList();

  /* get neibour list of atoms */
  cnfModifier.wrapAtomPrl(c0);
  cnfModifier.getNBL(c0, RCut);


#ifdef DEBUG
  // for (int i = 0; i < vacList.size(); ++i) {
  //   cout << vacList[i] << endl;
  //   for (int j = 0; j < c0.atoms[i].NBL.size(); ++j)
  //     cout << c0.atoms[i].NBL[j] << " ";
  //   cout << endl;
  // }
#endif
}



template<class ForwardIt, class T, class Compare>
ForwardIt mylower_bound(ForwardIt first, ForwardIt last, const T& value,\
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

  if(it != eventList.cend())
    cout << "prob: " << randVal << " count: " << distance(eventList.begin(), it)\
         << " event cprob: " << it->getcProb() << endl;
  return *it;
}

void KNHome::KMCSimulation() {
  KMCInit();
  step = 0;
  while (step < maxIter) {
    buildEventList();
    auto&& event = selectEvent();
    event.exeEvent(c0);
    ++step;
  }

}

void KNHome::buildEventList() {
  eventList.clear();
  double sum = 0.0;
  for (int i = 0; i < vacList.size(); ++i) {
    for (int j = 0; j < c0.atoms[i].NBL.size(); ++j) {
      KMCEvent event(std::make_pair(i, j));
      event.calRate(c0, temperature);
      sum += event.getRate();
      eventList.push_back(event);
    }
  }
  cout << "eventList size: " << eventList.size() << endl;

  double curr = 0.0;
  for (int i = 0; i < eventList.size(); ++i) {
    auto&& event = eventList[i];
    event.calProb(sum);
    curr += event.getProb();
    event.setcProb(curr);
  }

#ifdef DEBUG
  for (int i = 0; i < eventList.size(); ++i) {
    const auto& event = eventList[i];
    cout << setprecision(16) << i << " rate: " << event.getRate() \
         << " cumulative prob: " << event.getcProb() << endl;
  }
  cout << endl;
#endif
  
}