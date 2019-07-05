#include "gbCnf.h"
#include "KNHome.h"
#define KP 16
const string PBE="/Users/mingfei/work/pot_old/potpaw_PBE/elements/";

inline void prepINCAR(const string path) {
  string fnm = path + "/INCAR";
  ofstream ofs(fnm, std::ofstream::out);
  ofs << "NWRITE = 2\n";
  ofs << "                 \n";    
  ofs << "PREC   = Accurate\n";
  ofs << "ISYM   = 2       \n";    
  ofs << "NELM   = 240     \n";    
  ofs << "NELMIN = 4       \n";    
  ofs << "                 \n";    
  ofs << "NSW    = 10000   \n";    
  ofs << "IBRION = 2       \n";    
  ofs << "POTIM  = 0.5     \n";    
  ofs << "ISIF   = 2       \n";    
  ofs << "                 \n";    
  ofs << "ISMEAR = 1       \n";    
  ofs << "SIGMA  = 0.4     \n";    
  ofs << "                 \n";    
  ofs << "IALGO  = 48      \n";    
  ofs << "LREAL  = AUTO    \n";    
  ofs << "ENCUT  = 450.00  \n";    
  ofs << "ENAUG  = 600.00  \n";    
  ofs << "EDIFF  = 1e-7    \n";    
  ofs << "ISPIN  = 1       \n";    
  ofs << "                 \n";    
  ofs << "LWAVE  = .FALSE. \n";    
  ofs << "LCHARG = .TRUE.  \n";    
  ofs << "                 \n";    
  ofs << "NPAR   = 8       \n";    
}

inline void prepKPOINTS(const string path, const vector<int>& dupFac) {
  string fnm = path + "/KPOINTS";
  ofstream ofs(fnm, std::ofstream::out);
  ofs << "Automatic mesh\n";
  ofs << "0             \n";
  ofs << "Monkhorst-Pack\n";
  //ofs << KP/dupFac[X] << "   " << KP/dupFac[Y] << "   " << KP/dupFac[Z] << "\n";
  ofs << "1    1    1   \n";
  ofs << "0.   0.   0.  \n";
}  

inline void prepSUBMIT(const string path) {
  string fnm = path + "/submit.sh";
  ofstream ofs(fnm, std::ofstream::out);
  ofs << "/data/submit/unix/submit vasp ver=5.3.5 ncpu=48 spool_files=yes \
    queue=nahpc_matls_lg  proj=VASP  input_dir=`pwd` \
    jid=Al_job output_dir=`pwd`";
}

inline void prepPOTCAR(const string path, const set<string> species) {
  string mkPOT = "cat ";
  for (const auto& ele : species) {
    mkPOT += (PBE + ele + "/POTCAR ");
  }
  mkPOT += (" > " + path + "/POTCAR");
  const char *cmkPOT = mkPOT.c_str();

  const int tmp_err = std::system(cmkPOT);
  if (-1 == tmp_err) {
    cout << "Error making POTCAR in " << path << "\n";
  }
}

void KNHome::prepVASPFiles(const string path, const vector<int>& dupFac, 
                           const set<string>& species) {
  prepINCAR(path);
  prepKPOINTS(path, dupFac);
  prepPOTCAR(path, species);
  prepSUBMIT(path);
}
