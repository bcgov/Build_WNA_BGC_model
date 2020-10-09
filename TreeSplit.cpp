#include <Rcpp.h>
#include <stdio.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
vector<int> calcDepth(vector<int> nodes, vector<int> lchild, vector<int> rchild){
  int n = nodes.size();
  int currNode;
  int idx;
  vector<int> depth(n);
  vector<int>::iterator it;
  depth[0] = 0;
  for(int i = 1; i < n; i++){
    currNode = nodes[i];
    it = lower_bound(lchild.begin(),lchild.begin()+i,currNode);
    idx = it - lchild.begin();
    if(*it != currNode){
      it = lower_bound(rchild.begin(),rchild.begin()+i,currNode);
      idx = it - rchild.begin();
    }
    depth[i] = depth[idx]+1;
  }
  return(depth);
}

// [[Rcpp::export]]
IntegerVector findSplit(vector<int> node, vector<int> lchild, vector<int> rchild, 
                        vector<int> initNodes, int cutoff){
  long currNd;
  long currSz;
  int maxIdx;
  IntegerVector res;
  std::vector<int>::iterator itMax;
  std::map<long, long> Elt;
  Elt.clear();
  for(long i = 0; i < initNodes.size(); i++){
    Elt[initNodes[i]] += 1;
  }
  std::map<long,long>::iterator best;
  while(true){//need to always choose maximum value node
    Rcpp::checkUserInterrupt();
    itMax = max_element(initNodes.begin(),initNodes.end());
    maxIdx = itMax - initNodes.begin();
    //Rcout << "maxI " << maxIdx << "\n";
    currNd = *itMax;
    //Rcout << currNd << "\n";
    for(int j = (currNd-1); j > 0; j--){
      //Rcout << "j " << j << "\n";
      if(j == 1){ //so doesn't have infinite loop
        return(res);
      }
      if(lchild[j] == currNd || rchild[j] == currNd){
        initNodes[maxIdx] = node[j];
        Elt[node[j]] += 1;
        best = std::max_element(Elt.begin(),Elt.end(),[] 
          (const std::pair<long,long>& a, const std::pair<long,long>& b)->bool{ return a.second < b.second; });
        currSz = best->second;
        if(currSz > cutoff - 10){
          res.push_back(best->first);
        }
        break;
      }
    }
    if(currSz >= cutoff){
      return(res);
    }
  }
  return(res);
}
