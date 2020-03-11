
// modified from geeksforgeeks.com, originally contributed by Satish Srinivas


// We can use stl container list as a double
// ended queue to store the cache keys, with
// the descending time of reference from front
// to back and a set container to check presence
// of a key. But to fetch the address of the key
// in the list using find(), it takes O(N) time.
// This can be optimized by storing a reference
//   (iterator) to each key in a hash map.
#include "LRUCache.h"

// namespace LR {
LRUCache::LRUCache() : ct(0) {
}

// Declare the size
LRUCache::LRUCache(const int& cSizeIn)
  : ct(0), cSize(cSizeIn)  {
}

void LRUCache::setSize(const int& cSizeIn) {
  cSize = cSizeIn;
}
// Refers key x with in the LRU cache
void LRUCache::add(const pair<string, double>& x) {
  // not present in cache
  if (ma.find(x.first) == ma.end()) {
    // cache is full
    if (dq.size() == cSize) {
      // delete least recently used element
      string last = dq.back().first;

      // Pops the last elmeent
      dq.pop_back();

      // Erase the last
      ma.erase(last);
    }
  } else {// present in cache
    dq.erase(ma[x.first]);
  }

  // update reference
  dq.push_front(x);
  ma[x.first] = dq.begin();
}

void LRUCache::add(const pair<vector<int>, double>& x) {
  char xCStr[NB];
  int i = 0;
  for (const auto& digit : x.first)
    xCStr[i++] = '0' + digit;
  xCStr[i] = '\0';
  this->add(make_pair(xCStr, x.second));
}

bool LRUCache::check(const string& x) const {
  if (ma.find(x) != ma.end())
    return true;
  else
    return false;
}

bool LRUCache::check(const vector<int>& x) const {
  char xCStr[NB];
  int i = 0;
  for (const auto& digit : x)
    xCStr[i++] = '0' + digit;
  xCStr[i] = '\0';
  return check(xCStr);
}

double LRUCache::getBarrier(const string& x) {
  ++ct;
  return ma[x]->second;
}

double LRUCache::getBarrier(const vector<int>& x) {
  char xCStr[NB];
  int i = 0;
  for (const auto& digit : x)
    xCStr[i++] = '0' + digit;
  xCStr[i] = '\0';
  return getBarrier(xCStr);;
}

int LRUCache::getSize() const {
  return cSize;
}

// Function to display contents of cache
void LRUCache::display() const{

  // Iterate in the deque and print
  // all the elements in it
  for (auto it = dq.begin(); it != dq.end(); it++)
    cout << (*it).first << " " << (*it).second << endl;
}

long long LRUCache::getCt() const {
  return ct;
}
// } // end namespace LR