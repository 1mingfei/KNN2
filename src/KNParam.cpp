#include "gbCnf.h"
#include "KNHome.h"

void KNHome::initParam() { readParam(); }

/**************************************************
 * read parameters
 **************************************************/
void KNHome::readParam() {
  ifstream fid(sparams["parfile"], std::ifstream::in);
  if (!fid.is_open()) cerr << " error opening " << sparams["parfile"] << endl;
  vector<string> segs;
  string buff;
  while (getline(fid, buff)) {
    segs.clear();
    split(buff, " ", segs);
    if (!segs[0].compare("mode")) {
      sparams[segs[0]] = segs[1];
    } else if (!segs[0].compare("method")) {
      sparams[segs[0]] = segs[1];
    } else if (!segs[0].compare("initconfig")) {
      sparams[segs[0]] = segs[1];
    } else if (!segs[0].compare("PairFile")) {
      sparams[segs[0]] = segs[1];
    } else if (!segs[0].compare("ListFile")) {
      sparams[segs[0]] = segs[1];
    } else if (!segs[0].compare("format")) {
      sparams[segs[0]] = segs[1];
    } else if (!segs[0].compare("kerasModelBarrier")) {
      sparams[segs[0]] = segs[1];
    } else if (!segs[0].compare("kerasModelEDiff")) {
      sparams[segs[0]] = segs[1];
    } else if (!segs[0].compare("POT")) {
      sparams[segs[0]] = segs[1];
    } else if (!segs[0].compare("LC")) {
      dparams[segs[0]] = stod(segs[1]);
    } else if (!segs[0].compare("RCut")) {
      dparams[segs[0]] = stod(segs[1]);
    } else if (!segs[0].compare("temperature")) {
      dparams[segs[0]] = stod(segs[1]);
    } else if (!segs[0].compare("prefix")) {
      dparams[segs[0]] = stod(segs[1]);
    } else if (!segs[0].compare("factors")) {
      vector<int> tmpVec;
      for (unsigned int i = 1; i < segs.size(); ++i) {
        tmpVec.push_back(stoi(segs[i]));
      }
      viparams[segs[0]] = tmpVec;
    } else if (!segs[0].compare("elems")) {
      vector<string> tmpVec;
      for (unsigned int i = 1; i < segs.size(); ++i) {
        tmpVec.push_back(segs[i]);
      }
      vsparams[segs[0]] = tmpVec;
    } else if (!segs[0].compare("nums")) {
      vector<int> tmpVec;
      for (unsigned int i = 1; i < segs.size(); ++i) {
        tmpVec.push_back(stoi(segs[i]));
      }
      viparams[segs[0]] = tmpVec;
    } else if (!segs[0].compare("elemsEffectOffset")) {
      vector<double> tmpVec;
      for (unsigned int i = 1; i < segs.size(); ++i) {
        tmpVec.push_back(stod(segs[i]));
      }
      vdparams[segs[0]] = tmpVec;
    } else if (!segs[0].compare("NConfigs")) {
      iparams[segs[0]] = stoll(segs[1]);
    } else if (!segs[0].compare("low")) {
      iparams[segs[0]] = stoll(segs[1]);
    } else if (!segs[0].compare("high")) {
      iparams[segs[0]] = stoll(segs[1]);
    } else if (!segs[0].compare("NBarriers")) {
      iparams[segs[0]] = stoll(segs[1]);
    } else if (!segs[0].compare("maxIter")) {
      iparams[segs[0]] = stoll(segs[1]);
    } else if (!segs[0].compare("nTallyConf")) {
      iparams[segs[0]] = stoll(segs[1]);
    } else if (!segs[0].compare("nTallyOutput")) {
      iparams[segs[0]] = stoll(segs[1]);
    } else if (!segs[0].compare("randSeed")) {
      iparams[segs[0]] = stoll(segs[1]);
    } else if (!segs[0].compare("startingStep")) {
      iparams[segs[0]] = stoll(segs[1]);
    } else if (!segs[0].compare("startingTime")) {
      dparams[segs[0]] = stod(segs[1]);
    } else if (!segs[0].compare("startingEnergy")) {
      dparams[segs[0]] = stod(segs[1]);
    } else if (!segs[0].compare("EDiff")) {
      sparams[segs[0]] = segs[1];
    } else if (!segs[0].compare("ECutoff")) {
      dparams[segs[0]] = stod(segs[1]);
    } else if (!segs[0].compare("trapStepCriteria")) {
      iparams[segs[0]] = stoll(segs[1]);
    } else if (!segs[0].compare("oneTrapTimeCriteria")) {
      dparams[segs[0]] = stod(segs[1]);
    } else if (!segs[0].compare("LS_output_cfg_Criteria")) {
      dparams[segs[0]] = stod(segs[1]);
    } else if (!segs[0].compare("numClustersKept")) {
      iparams[segs[0]] = stoll(segs[1]);
    } else if (!segs[0].compare("solventAtomType")) {
      sparams[segs[0]] = segs[1];
    } else if (!segs[0].compare("numDataset")) {
      iparams[segs[0]] = stoll(segs[1]);
    } else if (!segs[0].compare("initNum")) {
      iparams[segs[0]] = stoll(segs[1]);
    } else if (!segs[0].compare("increment")) {
      iparams[segs[0]] = stoll(segs[1]);
    } else if (!segs[0].compare("finalNum")) {
      iparams[segs[0]] = stoll(segs[1]);
    } else if (!segs[0].compare("switchLSKMC")) {
      bparams[segs[0]] = (segs[1] == "false") ? false : true;
    } else if (!segs[0].compare("LRUSize")) {
      iparams[segs[0]] = stoll(segs[1]);
    } else if (!segs[0].compare("solventBoudCriteria")) {
      iparams[segs[0]] = stoll(segs[1]);
    } else if (!segs[0].compare("smallestClusterCriteria")) {
      iparams[segs[0]] = stoll(segs[1]);
    }
  }
  fid.close();
}

/**************************************************
 * parse arguments
 **************************************************/
void KNHome::parseArgs(int argc, char* argv[]) {
  for (int i = 0; i < argc; ++i) {
    if (!strcmp(argv[i], "--p") || !strcmp(argv[i], "-p"))
      sparams["parfile"] = string(argv[++i]);
    if (!strcmp(argv[i], "--i") || !strcmp(argv[i], "-i"))
      sparams["dmpfile"] = string(argv[++i]);
    if (!strcmp(argv[i], "--f") || !strcmp(argv[i], "-f"))
      sparams["potfile"] = string(argv[++i]);
  }
}
