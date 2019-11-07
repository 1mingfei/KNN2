#include "KMCEvent.h"

#define KB 8.6173303e-5

KMCEvent::KMCEvent(const pair<int, int>& inPair)
: jumpPair(inPair) {}

KMCEvent::KMCEvent() {
}

KMCEvent::~KMCEvent() {}

void KMCEvent::getRate(const Config& cnf) {

}

void KMCEvent::calRate(const double& deltaE, const double& T) {
  rate = exp(deltaE / KB / T);
}

void KMCEvent::exeEvent(Config& cnf) {

}