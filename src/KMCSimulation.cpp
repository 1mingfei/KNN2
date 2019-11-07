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
  c0 = cnfModifier.readCfg(fname);

  /* get initial vacancy position in atomList */
  getVacList();

  /* get neibour list of atoms */
  cnfModifier.wrapAtomPrl(c0);
  cnfModifier.getNBL(c0, RCut);


#ifdef DEBUG
  for (int i = 0; i < vacList.size(); ++i) {
    cout << vacList[i] << endl;
    for (int j = 0; j < c0.atoms[i].NBL.size(); ++j)
      cout << c0.atoms[i].NBL[j] << " ";
    cout << endl;
  }
#endif


}

void KNHome::KMCSimulation() {
  KMCInit();
  int step = 0;
  while (step < maxIter) {
    buildEventList();
  }

}

void KNHome::buildEventList() {

  for (int i = 0; i < vacList.size(); ++i) {
    for (int j = 0; j < c0.atoms[i].NBL.size(); ++j) {
      KMCEvent event(std::make_pair(i, j));
      eventList.push_back(event);
    }
  }
  
}