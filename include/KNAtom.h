#ifndef _KNATOM_
#define _KNATOM_

using std::vector;
using std::string;

class KNAtom {
private:
public:
  int id;
  string tp;
  double pst[3], prl[3];
  vector<int> NBL;

  KNAtom() : id(0), tp("H") {};
  KNAtom(int n) : id(n), tp("H") {};
  KNAtom(int _id, double x, double y, double z) : id(_id), tp("H") {
    prl[0] = x, prl[1] = y, prl[2] = z;
  };

  KNAtom(int _id, string _tp, double x, double y, double z):
    id(_id), tp(_tp) {
    prl[0] = x, prl[1] = y, prl[2] = z;
  }

  ~KNAtom(){};

  bool operator<(const KNAtom &b) const { return tp < b.tp; }
};

#endif
