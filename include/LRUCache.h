#ifndef _LRU_CACHE_H_
#define _LRU_CACHE_H_
// modified from geeksforgeeks.com, originally contributed by Satish Srinivas


// We can use stl container list as a double
// ended queue to store the cache keys, with
// the descending time of reference from front
// to back and a set container to check presence
// of a key. But to fetch the address of the key
// in the list using find(), it takes O(N) time.
// This can be optimized by storing a reference
//   (iterator) to each key in a hash map.

#include "KNHome.h"
#define NB 28
#define KEY_SIZE 27
class gbCnf;

namespace std
{
template<typename T>
struct hash<vector<T> > {
  typedef vector<T> argument_type;
  typedef size_t result_type;
  size_t N = KEY_SIZE;

  result_type operator()(const argument_type& a) const {
    result_type h = a[0];
    for (result_type i = 1; i < N; ++i) {
      h <<= 2;
      h += a[i];
    }
    return h;
  }
};
}

// namespace LR {

using std::list;

class LRUCache {
  long long ct;
  // store keys of cache
  list<pair<string, double>> dq;

  // store references of key in cache
  unordered_map<string, list<pair<string, double>>::iterator> m;
  int cSize; // maximum capacity of cache

public:
  LRUCache();
  LRUCache(const int&);
  void setSize(const int&);
  void add(const pair<vector<int>, double>&);
  void add(const pair<string, double>&);
  bool check(const vector<int>&) const;
  bool check(const string&) const;
  double getBarrier(const vector<int>&);
  double getBarrier(const string&);
  int getSize() const;
  void display();
  long long getCt() const;
};
// } // end namespace LR
#endif