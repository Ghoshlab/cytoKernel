#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export(.getMatchIndices2)]]
IntegerVector getMatchIndices2(CharacterVector x, CharacterVector match) {
  std::vector<int> idx;
  for (int i = 0; i < x.size(); ++i) {
    if (x[i] == match[0]) {
      idx.push_back(i + 1); // add 1 to match R's indexing
    }
  }
  return wrap(idx);
}





