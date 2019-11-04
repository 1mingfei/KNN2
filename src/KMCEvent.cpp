#include "KMCEvent.h"
#include "KNHome.h"

KNHome::KMCEvent::KMCEvent(KNHome& x)
 :kn(x),
  dparams(x.dparams),
  iparams(x.iparams), 
  sparams(x.sparams), 
  vsparams(x.vsparams), 
  viparams(x.viparams) 
{}
