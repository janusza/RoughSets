#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List compute_indiscernibility(List input, CharacterVector attr_val) {
  std::vector<IntegerVector> result;
  int input_size = input.size();
  for (int i=0; i<input_size; ++i) {
    IntegerVector ind_class = input[i];
    std::map<String const, std::vector<int> > split;
    int ind_class_size = ind_class.size();
    for (int j=0; j<ind_class_size; ++j) {
      split[attr_val[ind_class[j]-1]].push_back(ind_class[j]);  
    }
    for (std::map<String const, std::vector<int> >::iterator it = split.begin(); it != split.end(); ++it) {
      result.push_back(wrap(it->second));
    }
  } 
  return wrap(result);
}
