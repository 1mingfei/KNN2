#include "KMCEvent.h"

#define KB 8.6173303e-5

KMCEvent::KMCEvent(const pair<int, int>& inPair)
: jumpPair(inPair) {}

KMCEvent::KMCEvent() {
}

KMCEvent::~KMCEvent() {}

double KMCEvent::getRate() const {
  return rate;
}

double KMCEvent::getProb() const {
  return prob;
}

double KMCEvent::getcProb() const {
  return cProb;
}

void KMCEvent::calRate(const Config& cnf, const double& T) {

  // need to be changed
  double deltaE = (double) rand() / (RAND_MAX);
  rate = exp(- deltaE / KB / T);
#ifdef DEBUG
  // cout << "activation barrier: " << deltaE << " rate: " << rate << endl;
#endif
}

void KMCEvent::calProb(const double& sum) {
  prob = rate / sum;
}

void KMCEvent::setcProb(const double& curr) {
  cProb = curr;
}


void KMCEvent::exeEvent(Config& cnf) {

}