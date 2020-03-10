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

class gbCnf;

// namespace LR {

using std::list;

class LRUCache {
  long long ct;
  // store keys of cache
  list<pair<string, double>> dq;

  // store references of key in cache
  unordered_map<string, list<pair<string, double>>::iterator> ma;
  int cSize; // maximum capacity of cache

public:
  LRUCache();
  LRUCache(const int&);
  void setSize(const int&);
  void add(const pair<vector<int>, double>&);
  void add(const pair<string, double>&);
  bool check(const string&) const;
  bool check(const vector<int>&) const;
  double getBarrier(const string&);
  double getBarrier(const vector<int>&);
  int getSize() const;
  void display() const;
  long long getCt() const;
};
// } // end namespace LR
#endif