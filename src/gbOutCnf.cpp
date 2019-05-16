/*
 * @Author: chaomy
 * @Date:   2018-07-07 16:58:27
 * @Last Modified by:  1mingfei 
 * @Last Modified time: 2019-5-16
 */

#include "gbCnf.h"

double KNHome::gbCnf::findMass(string x) {
  /*this "it" is indeed a iterator*/
  auto it = std::find(element.begin(), element.end(), x);
  int index = std::distance(element.begin(), it);
  return mass[index];
}
/**************************************************
 * write lammps data files
 **************************************************/
/*
  #lmp data config
  912 atoms
  3 atom types
  0.000000  31.772283 xlo xhi
  0.000000  200.000000 ylo yhi
  0.000000  11.071269 zlo zhi
  0.000000  0.000000  0.000000 xy xz yz
  Atoms

*/
void KNHome::gbCnf::writeLmpData(Config& c, string fnm = "out.lmp.init") {
  ofstream ofs(fnm, std::ofstream::out);
  ofs << "# lmp data config " << endl;
  ofs << (int)c.natoms << " atoms" << endl;
  //ofs << (int)c.ntypes << " atom types " << endl;
  ofs << 4 << " atom types " << endl;
  ofs << c.cell[0] << " " << c.cell[3] << " xlo xhi" << endl;
  ofs << c.cell[1] << " " << c.cell[4] << " ylo yhi" << endl;
  ofs << c.cell[2] << " " << c.cell[5] << " zlo zhi" << endl;
  ofs << c.cell[6] << " " << c.cell[7] << " " << c.cell[8] << " xy xz yz"
      << endl;
  ofs << "Atoms" << endl << endl;
  for (int i = 0; i < c.natoms; ++i) {
    auto&& a = c.atoms[i];
    ofs << a.id + 1 << " " << a.tp << " " << a.pst[0] << " " << a.pst[1]
        << " " << a.pst[2] << endl;
  }
}

/* write cfg data files*/
/*
Number of particles = 16200
A = 4.37576470588235 Angstrom (basic length-scale)
H0(1,1) = 127.5 A
H0(1,2) = 0 A
H0(1,3) = 0 A
H0(2,1) = 0 A
H0(2,2) = 119.501132067411 A
H0(2,3) = 0 A
H0(3,1) = 0 A
H0(3,2) = 0 A
H0(3,3) = 3 A
.NO_VELOCITY.
entry_count = 9
auxiliary[0] = kine [reduced unit]
auxiliary[1] = pote [reduced unit]
auxiliary[2] = s11 [reduced unit]
auxiliary[3] = s22 [reduced unit]
auxiliary[4] = s12 [reduced unit]
auxiliary[5] = hydro [reduced unit]
1.000000
Ar
0.0016667 0.00616 0.5 0 -2 -1.9431e-13 -2.9917e-13 1.6811e-13 -2.4674e-13 

      buffData[0] = meanX;
      buffData[1] = stdX;
      buffData[2] = meanY;
      buffData[3] = stdY;
      buffData[4] = meanZ;
      buffData[5] = stdZ;
      buffData[6] = meanDist;
      buffData[7] = stdDist;
      buffData[8] = meanType;
      buffData[9] = stdTp;

 */
void KNHome::gbCnf::writeCfgData(const Config& c, 
                                 string fnm = "out.cfg") {
  ofstream ofs(fnm, std::ofstream::out);
  ofs << "Number of particles =  " << c.natoms << endl;
  ofs << "A = 1.0 Angstrom (basic length-scale)" << endl;
  ofs << "H0(1,1) = " << c.bvx[X] << " A" << endl;
  ofs << "H0(1,2) = " << c.bvx[Y] << " A" << endl;
  ofs << "H0(1,3) = " << c.bvx[Z] << " A" << endl;
  ofs << "H0(2,1) = " << c.bvy[X] << " A" << endl;
  ofs << "H0(2,2) = " << c.bvy[Y] << " A" << endl;
  ofs << "H0(2,3) = " << c.bvy[Z] << " A" << endl;
  ofs << "H0(3,1) = " << c.bvz[X] << " A" << endl;
  ofs << "H0(3,2) = " << c.bvz[Y] << " A" << endl;
  ofs << "H0(3,3) = " << c.bvz[Z] << " A" << endl;
  ofs << ".NO_VELOCITY." << endl;
  ofs << "entry_count = 3" << endl;
  double mass0 = findMass(c.atoms[0].tp);
  //double mass1 = findMass(elems[1]);
  for (int i = 0; i < c.atoms.size(); ++i) {
    auto&& a = c.atoms[i];
    ofs << mass0 << "\n" << a.tp << "\n";
    ofs << a.prl[X] << " " << a.prl[Y] << " " << a.prl[Z] << "\n";
  }
}
/*
 #Comment
1.000000
     16.18400000        0.00000000        0.00000000
      0.00000000       16.18400000        0.00000000
      0.00000000        0.00000000       16.18400000
      A      B
   240      13
Cartesian
      0.00000000        0.00000000        0.00000000
     12.13800000       14.16100000       14.16100000
 */
/* this output vac as X to visualize it */
void KNHome::gbCnf::writePOSCARVis(Config& c, string fnm = "POSCAR", \
                                   string comment = "") {
  ofstream ofs(fnm, std::ofstream::out);
  ofs << "# " << comment << "\n";
  ofs << "1.00000\n";
  for (int i = 0; i < DIM; ++i) {
    ofs << c.bvx[i] << " ";
  }
  ofs << "\n";
  for (int i = 0; i < DIM; ++i) {
    ofs << c.bvy[i] << " ";
  }
  ofs << "\n"; 
  for (int i = 0; i < DIM; ++i) {
    ofs << c.bvz[i] << " ";
  }
  ofs << "\n";
  std::map<string, int> names;
  for (const auto & atm : c.atoms) {
    if (names.find(atm.tp) == names.end()) {
      names[atm.tp] = 1;
    } else {
      ++names[atm.tp];
    }
  }
  for (auto i = names.begin(); i != names.end(); ++i) {
    ofs << i->first << " "; 
  }
  ofs << "\n";
  for (auto i = names.begin(); i != names.end(); ++i) {
    ofs << i->second << " "; 
  }
  //std::sort(c.atoms.begin(), c.atoms.end());
  ofs << "\nDirect\n";
  for (unsigned int i = 0 ; i < c.atoms.size(); ++i) {
    for (int j = 0; j < DIM; ++j) {
      ofs << c.atoms[i].prl[j] << " "; 
    }
    ofs << "\n";
  }
}

/* this output is for vasp calculation (ignore vac) */
void KNHome::gbCnf::writePOSCAR(Config& c, string fnm = "POSCAR") {
  ofstream ofs(fnm, std::ofstream::out);
  ofs << "#comment\n1.00000\n";
  for (int i = 0; i < DIM; ++i) {
    ofs << c.bvx[i] << " ";
  }
  ofs << "\n";
  for (int i = 0; i < DIM; ++i) {
    ofs << c.bvy[i] << " ";
  }
  ofs << "\n"; 
  for (int i = 0; i < DIM; ++i) {
    ofs << c.bvz[i] << " ";
  }
  ofs << "\n";
  std::map<string, int> names;
  for (const auto & atm : c.atoms) {
    if (names.find(atm.tp) == names.end()) {
      names[atm.tp] = 1;
    } else {
      ++names[atm.tp];
    }
  }
  for (auto i = names.begin(); i != names.end(); ++i) {
    if (i->first != "X") {
      ofs << i->first << " ";
    }
  }
  ofs << "\n";
  for (auto i = names.begin(); i != names.end(); ++i) {
    if (i->first != "X") {
      ofs << i->second << " ";
    }
  }
  //std::sort(c.atoms.begin(), c.atoms.end());
  ofs << "\nDirect\n";
  for (unsigned int i = 0 ; i < c.atoms.size(); ++i) {
    if (c.atoms[i].tp != "X") {
      for (int j = 0; j < DIM; ++j) {
        ofs << c.atoms[i].prl[j] << " "; 
      }
      ofs << "\n";
    }
  }
}
