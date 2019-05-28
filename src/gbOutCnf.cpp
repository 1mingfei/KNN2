/*
 * @Author: chaomy
 * @Date:   2018-07-07 16:58:27
 * @Last Modified by:  1mingfei 
 * @Last Modified time: 2019-5-26
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

/* write cfg data files */
void KNHome::gbCnf::writeCfgData(const Config& c, 
                                 string fnm = "out.cfg") {
  ofstream ofs(fnm, std::ofstream::out);
  ofs << "Number of particles = " << c.natoms << endl;
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
  double mass = findMass(c.atoms[0].tp);
  string prevType = c.atoms[0].tp;
  for (int i = 0; i < c.atoms.size(); ++i) {
    auto&& a = c.atoms[i];
    if (a.tp != prevType) {
      mass = findMass(a.tp);
    }
    ofs << mass << "\n" << a.tp << "\n";
    ofs << a.prl[X] << " " << a.prl[Y] << " " << a.prl[Z] << "\n";
  }
}

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
