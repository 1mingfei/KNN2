#ifndef _GB_DEF_
#define _GB_DEF_

#include "KNAtom.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::pair;

#define DIM 3
#define MAXLEN 1024
#define PI 3.14159265359
#define INVPI 0.31830988618
#define SQRT3 1.73205080757
#define KB 8.6173303e-5

enum { XX = 0, YY = 1, ZZ = 2, XY = 3, YZ = 4, ZX = 5 };
enum { X = 0, Y = 1, Z = 2 };

class Config {
public:
  int natoms, ntypes;
  double engy;

  vector<double> cell;    // lox, loy, loz, hix, hiy, hiz, xy xz yz
  vector<double> length;  // length of three edges
  vector<double> bvx, tvx, bvy, tvy, bvz, tvz;
  vector<KNAtom> atoms;
  vector<int> vacList;

  bool operator<(const Config &b) const { return this->engy < b.engy; }
  Config()
      : natoms(0),
        ntypes(0),
        engy(0.0),
        cell(9),
        length(3),
        bvx(3),
        tvx(3),
        bvy(3),
        tvy(3),
        bvz(3),
        tvz(3){};
  ~Config(){};
friend class KNHome;
};

#endif
