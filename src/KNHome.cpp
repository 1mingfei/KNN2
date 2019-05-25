/*
 * Author: 1mingfei 
 * Date:   2019-04-14
 * Purpose: inline functions for KNHome
 * self-explained
 */

#include "gbCnf.h"
#include "KNHome.h"

KNHome::KNHome(int argc, char* argv[]) {
  MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  sparams["mode"] = "generate";
  dparams["LC"] = 4.046;
  dparams["Rcut"] = 6.0;

  parseArgs(argc, argv);
  srand(time(NULL) + me);
  initParam();

  gbCnf cnfModifier(*this);
  
  if (sparams["mode"] == "generate") {
    createPreNEB();
  } else if (sparams["mode"] == "encode") {
    //KNEncode();
  }

}

KNHome::~KNHome() {}
