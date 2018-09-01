#ifndef _disjoint_dense_int_sets_
#define _disjoint_dense_int_sets_

#include <vector>
#include <cassert>

class DisjointDenseIntSet {
 private:
  std::vector<unsigned int> rank;
  std::vector<unsigned int> parent;
  void checkSize(unsigned int newN){
    if(newN<rank.size())
      return;

    const auto old_size = rank.size();
    rank.resize(newN+1,0);
    parent.resize(newN+1);
    for(auto i=old_size;i<newN+1;i++)
      parent[i] = i;
  }
 public:
  DisjointDenseIntSet(){}
  DisjointDenseIntSet(unsigned int N){
    rank.resize(N,0);
    parent.resize(N);
    for(unsigned int i=0;i<N;i++)
      parent[i] = i;
  }
  void makeSet(unsigned int n){
    checkSize(n);
    parent[n] = n;
    rank[n]   = 0;
  }
  unsigned int maxElement() const {
    return rank.size()-1;
  }
  unsigned int findSet(unsigned int n){
    if(parent[n]==n)
      return n;
    else
      return parent[n] = findSet(parent[n]);
  }
  void unionSet(unsigned int a, unsigned int b){
    auto roota = findSet(a);
    auto rootb = findSet(b);
    if(roota==rootb)
      return;

    if(rank[roota]<rank[rootb]){
      parent[roota] = rootb;
    } else if(rank[roota]>rank[rootb]) {
      parent[rootb] = roota;
    } else {
      parent[rootb] = roota;
      rank[roota]++;
    }
  }
  void mergeAintoB(unsigned int a, unsigned int b){
    checkSize(a);
    checkSize(b);
    parent[a] = b;
  }
  bool sameSet(unsigned int a, unsigned int b){
    return findSet(a)==findSet(b);
  }
};

#endif