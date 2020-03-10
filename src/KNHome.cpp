#include "gbCnf.h"
#include "KNHome.h"
#include "LSKMC.h"
#include "LRUCache.h"

KNHome::KNHome(int argc, char* argv[]) {
  MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  sparams["mode"] = "generate";
  sparams["initconfig"] = "initconfig.cfg";
  dparams["LC"] = 4.046;
  // dparams["Rcut"] = 6.0; //encode
  dparams["Rcut"] = 3.0; //kmc
  iparams["randSeed"] = 1234567; //kmc

  gbCnf cnfModifier(me, nProcs);

  if (!strcmp(argv[1], "POSCAR") && (me == 0)) {

    Config cfg = std::move(cnfModifier.readPOSCAR("POSCAR"));

    string mkBaseDir = "mv POSCAR POSCAR.ori";
    const char *cmkBaseDir = mkBaseDir.c_str();
    const int mv_err = std::system(cmkBaseDir);
    if (-1 == mv_err) {
      cout << "Error moving original POSCAR\n";
      exit(1);
    }

    cnfModifier.perturb(cfg);
    map<string, int> elemName = cnfModifier.writePOSCAR(cfg, "POSCAR");
    exit(0);
  }

  parseArgs(argc, argv);
  initParam();
  if (me == 0) {
    LRUSize = iparams["LRUSize"];
    cout << "LRUSize set to : " << LRUSize << endl;
    cout << "LRUCache will " << (LRUSize ? "" : "not ") << "be used." << endl;
  }
  lru = new LRUCache(LRUSize);

  if ((sparams["mode"]) == "generate") {
    srand(iparams["randSeed"] + me);
    createPreNEB(cnfModifier);
  } else if (sparams["mode"] == "encode") {
    KNEncode(cnfModifier);
  } else if (sparams["mode"] == "BondCount") {
    KNBondCount(cnfModifier);
  } else if (sparams["mode"] == "kmc") {
    KMCSimulation(cnfModifier);
  } else if (sparams["mode"] == "lskmc_test") {
    LSKMCOneRun(cnfModifier);
  } else if (sparams["mode"] == "lskmc") {
    LSKMCSimulation(cnfModifier);
  } else if (sparams["mode"] == "clusterCount") {
    loopConfig(cnfModifier);
  } else if (sparams["mode"] == "test") {
    // testK2P();
    /* test encoding */
    // gbCnf cnfModifier(*this);
    // Config c1 = cnfModifier.readCfg("in.cfg");
    // cnfModifier.writePOSCARVis(c1, "POSCAR", "");
    // cnfModifier.writeCfgData(c1, "test_out.cfg");
    /* test LSKMC */
    // LSKMCSimulation(cnfModifier);

    /* test mat transfer */
    // LS::LSKMC::test_vvd2mat();
    /* test embedding */
    buildEmbedding();


  }
}

KNHome::~KNHome() {
  delete lru;
}
