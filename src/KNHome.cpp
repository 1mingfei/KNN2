/*
 * Author: 1mingfei 
 * Date:   2019-05-27
 * Purpose: functions for KNHome
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
  
  if ((sparams["mode"]) == "generate") {
    createPreNEB();
  } else if (sparams["mode"] == "encode") {
    KNEncode();
  } else if (sparams["mode"] == "BondCount") {
    KNBondCount();
  // } else if (sparams["mode"] == "kmc") {
  //   KMCSimulation();
  } else if (sparams["mode"] == "test") {
    /* test encoding */
    /*
    gbCnf cnfModifier(*this);
    Config c1 = cnfModifier.readCfg("in.cfg");
    cnfModifier.writePOSCARVis(c1, "POSCAR", "");
    cnfModifier.writeCfgData(c1, "test_out.cfg");
    */
    /* test uniform type */
  }


}

KNHome::~KNHome() {}
