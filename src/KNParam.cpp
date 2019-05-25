/*
 * @Author: chaomy
 * @Date:   2017-12-31 16:04:15
 * @Last Modified by:  1mingfei
 * @Last Modified time: 2019-03-25
 */

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
    } else if (!segs[0].compare("LC")) {
      dparams[segs[0]] = stod(segs[1]);
    } else if (!segs[0].compare("RCut")) {
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
    } else if (!segs[0].compare("NConfigs")) {
      iparams[segs[0]] = stoi(segs[1]);
    } else if (!segs[0].compare("NBarriers")) {
      iparams[segs[0]] = stoi(segs[1]);
    }

  }
  fid.close();
}

/**************************************************
 * parse arguments
 **************************************************/
void KNHome::parseArgs(int argc, char* argv[]) {
  for (int i = 0; i < argc; i++) {
    if (!strcmp(argv[i], "--p") || !strcmp(argv[i], "-p"))
      sparams["parfile"] = string(argv[++i]);
    if (!strcmp(argv[i], "--i") || !strcmp(argv[i], "-i"))
      sparams["dmpfile"] = string(argv[++i]);
    if (!strcmp(argv[i], "--f") || !strcmp(argv[i], "-f"))
      sparams["potfile"] = string(argv[++i]);
  }
}