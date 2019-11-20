#include "gbCnf.h"
#include "KNHome.h"
#define KP 9
const string PBE="/Users/mingfei/work/pot_old/potpaw_PBE/elements/";
//const string PBE="/home/mingfei/Work/pot_old/potpaw_PBE/elements/";

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
  ofs << "NPAR   = 4       \n";    
}

inline void prepKPOINTS(const string path, const vector<int>& dupFac) {
  string fnm = path + "/KPOINTS";
  ofstream ofs(fnm, std::ofstream::out);
  ofs << "Automatic mesh\n";
  ofs << "0             \n";
  ofs << "Monkhorst-Pack\n";
  ofs << KP/dupFac[X] << "   " << KP/dupFac[Y] << "   " << KP/dupFac[Z] << "\n";
  //ofs << "1    1    1   \n";
  ofs << "0.   0.   0.  \n";
}  

inline void prepSUBMIT(const string path) {
  string fnm = path + "/submit.sh";
  ofstream ofs(fnm, std::ofstream::out);
  ofs << "/data/submit/unix/submit vasp ver=5.3.5 ncpu=48 spool_files=yes \
    queue=nahpc_matls_lg  proj=VASP  input_dir=`pwd` \
    jid=Al_job output_dir=`pwd`";
}

inline void prepSUBMITCORI(const string path) {
  string fnm = path + "/submit.cori";
  ofstream ofs(fnm, std::ofstream::out);
  ofs << "#!/bin/bash\n";
  ofs << "#SBATCH -N 1\n";
  ofs << "#SBATCH -C knl\n";
  ofs << "#SBATCH -q regular\n";
  ofs << "#SBATCH -J start\n";
  ofs << "#SBATCH --mail-user=mingfei@umich.edu\n";
  ofs << "#SBATCH --mail-type=ALL\n";
  ofs << "#SBATCH -t 48:00:00\n";
  ofs << "#SBATCH -L SCRATCH\n";
  ofs << "\n";
  ofs << "#OpenMP settings:\n";
  ofs << "export OMP_NUM_THREADS=4\n";
  ofs << "export OMP_PLACES=threads\n";
  ofs << "export OMP_PROC_BIND=spread\n";
  ofs << "\n";
  ofs << "module load vasp/5.4.4-knl\n";
  ofs << "srun -n 64 -c 4 --cpu_bind=cores vasp_std\n";
  ofs << "rm CHG* WAVE*\n";
}  

inline void prepSUBMITGL(const string path) {
  string fnm = path + "/submit.gl";
  ofstream ofs(fnm, std::ofstream::out);

  ofs << "#!/bin/bash\n";
  ofs << "#SBATCH -N 1\n";
  ofs << "#SBATCH --ntasks-per-node=36\n";
  ofs << "#SBATCH --mail-user=mingfei@umich.edu\n";
  ofs << "#SBATCH --mail-type=ALL\n";
  ofs << "#SBATCH -t 96:00:00\n";
  ofs << "#SBATCH --job-name=GOALI\n";
  ofs << "#SBATCH --account=qiliang\n";
  ofs << "#SBATCH --partition=standard\n";
  ofs << "\n";
  ofs << "module load RestrictedLicense\n";
  ofs << "module load vasp/5.4.4.18Apr17.p1\n";
  ofs << "srun  vasp\n";
  ofs << "rm CHG* WAVE*\n";
}

inline void prepSUBMITSTAMPEDE2(const string path) {
  string fnm = path + "/submit.stampede2";
  ofstream ofs(fnm, std::ofstream::out);
  ofs << "#!/bin/bash\n";
  ofs << "#SBATCH -J vasp\n";
  ofs << "#SBATCH -o vasp.%j.out\n";
  ofs << "#SBATCH -e vasp.%j.err\n";
  ofs << "#SBATCH -n 64\n";
  ofs << "#SBATCH -N 1\n";
  ofs << "#SBATCH -p normal\n";
  ofs << "#SBATCH -t 32:00:00\n";
  ofs << "#SBATCH -A  TG-MSS160003\n";
  ofs << "\n";
  ofs << "#TG-MSS160003 TG-DMR190035\n";
  ofs << "\n";
  ofs << "module load vasp/5.4.4\n";
  ofs << "ibrun vasp_std > vasp_test.out\n";
  ofs << "rm CHG* WAVE*\n";
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

void KNHome::prepVASPFiles(const string path, \
                           const vector<int>& dupFac, \
                           const set<string>& species) {
  prepINCAR(path);
  prepKPOINTS(path, dupFac);
  prepPOTCAR(path, species);
  prepSUBMIT(path);
  prepSUBMITCORI(path);
  prepSUBMITGL(path);
  prepSUBMITSTAMPEDE2(path);
}
