#include <Rcpp.h>
using namespace Rcpp;

template<typename T>
inline bool filter(T const & x){
  return x.size() <= 1;
}

// [[Rcpp::export]]
List compute_indiscernibility(List input, CharacterVector attr_val, CharacterVector unique_attr_val) {
  std::map<String, int> numbers;
  for (int i=0; i<unique_attr_val.size(); ++i) {
    numbers[unique_attr_val[i]] = i;
  }
  int numbers_size = numbers.size();
  int input_size = input.size();
  std::vector<std::vector<int> > result(input_size*numbers_size);
  for (int i=0; i<input_size; ++i) {
    IntegerVector ind_class = input[i];
    int ind_class_size = ind_class.size();
    for (int j=0; j<ind_class_size; ++j) {
      result[i*numbers_size+numbers[attr_val[ind_class[j]-1]]].push_back(ind_class[j]);
    }
  }
  result.erase(std::remove_if(result.begin(), result.end(), filter<std::vector<int> >), result.end());
  return wrap(result);
}

// [[Rcpp::export]]
List compute_chaos(List input, CharacterVector dec_val, CharacterVector unique_dec_val) {

  std::map<String, int> indexes;
  for (int i=0; i<unique_dec_val.size(); ++i) {
    indexes[unique_dec_val[i]] = i;
  }
  int input_size = input.size();
  std::vector<std::vector<int> > result(input_size);
  int unique_dec_val_size = unique_dec_val.size();
  std::vector<int> classCounts(unique_dec_val_size);

  for (int i=0; i<input_size; ++i) {
    for (int j=0; j<unique_dec_val_size; ++j) {
      classCounts[j] = 0;
    }
    IntegerVector ind_class = input[i];
    int ind_class_size = ind_class.size();
    for (int j=0; j<ind_class_size; ++j) {
      ++classCounts[indexes[dec_val[ind_class[j]-1]]];
    }
    result[i] = classCounts;
  }
  return wrap(result);
}

// // [[Rcpp::export]]
// List compute_indiscernibility_and_chaos(List input, CharacterVector attr_val, CharacterVector unique_attr_val, CharacterVector dec_val, CharacterVector unique_dec_val) {
//
//   std::map<String, int> numbers;
//   for (int i=0; i<unique_attr_val.size(); ++i) {
//     numbers[unique_attr_val[i]] = i;
//   }
//   int numbers_size = numbers.size();
//   int input_size = input.size();
//
//   std::vector<std::vector<int> > new_ind_rel(input_size*numbers_size);
//   for (int i=0; i<input_size; ++i) {
//     IntegerVector ind_class = input[i];
//     int ind_class_size = ind_class.size();
//     for (int j=0; j<ind_class_size; ++j) {
//       new_ind_rel[i*numbers_size+numbers[attr_val[ind_class[j]-1]]].push_back(ind_class[j]);
//     }
//   }
//   new_ind_rel.erase(std::remove_if(new_ind_rel.begin(), new_ind_rel.end(), filter<std::vector<int> >), new_ind_rel.end());
//
//   std::map<String, int> indexes;
//   for (int i=0; i<unique_dec_val.size(); ++i) {
//     indexes[unique_dec_val[i]] = i;
//   }
//   int new_ind_rel_size = new_ind_rel.size();
//
//   std::vector<std::vector<int> > myresult;
//   if (new_ind_rel_size > 0) {
//     myresult.resize(new_ind_rel_size + 1);
//     int unique_dec_val_size = unique_dec_val.size();
//     std::vector<int> classCounts(unique_dec_val_size);
//     for (int i=0; i<new_ind_rel_size; ++i) {
//       for (int j=0; j<unique_dec_val_size; ++j) {
//         classCounts[j] = 0;
//       }
//       std::vector<int> ind_class = new_ind_rel[i];
//       int ind_class_size = ind_class.size();
//       myresult[0].push_back(ind_class_size);
//       for (int j=0; j<ind_class_size; ++j) {
//         ++classCounts[indexes[dec_val[ind_class[j]-1]]];
//       }
//       myresult[i+1] = classCounts;
//     }
//   } else {
//     myresult.resize(2);
//     myresult[0].push_back(10);
//     myresult[1].push_back(10);
//     myresult[1].push_back(0);
//   }
//   return wrap(myresult);
// }

//// [[Rcpp::export]]
// List compute_indiscernibility_and_chaos2(List input, CharacterVector attr_val, CharacterVector unique_attr_val, CharacterVector dec_val, CharacterVector unique_dec_val) {
//
//   std::map<String, int> numbers;
//   int unique_attr_val_size = unique_attr_val.size();
//   for (int i=0; i<unique_attr_val_size; ++i) {
//     numbers[unique_attr_val[i]] = i;
//   }
//
//   std::map<String, int> indexes;
//   int unique_dec_val_size = unique_dec_val.size();
//   for (int i=0; i<unique_dec_val_size; ++i) {
//     indexes[unique_dec_val[i]] = i;
//   }
//
//   std::vector<int> tmp_class_counts(unique_dec_val_size, 0);
//
//   int input_size = input.size();
//   std::vector<std::vector<int> > class_counts_list(input_size*unique_attr_val_size, tmp_class_counts);
//   std::vector<int> final_ind_classes_size(input_size*unique_attr_val_size, 0);
//   int tmp_number = 0;
//   for (int i=0; i<input_size; ++i) {
//     std::vector<int> ind_class = input[i];
//     int ind_class_size = ind_class.size();
//
//     for (int j=0; j<ind_class_size; ++j) {
//       tmp_number = i*unique_attr_val_size + numbers[attr_val[ind_class[j]-1]];
//       ++final_ind_classes_size[tmp_number];
//       ++class_counts_list[tmp_number][indexes[dec_val[ind_class[j]-1]]];
//     }
//   }
//
//   class_counts_list.insert(class_counts_list.begin(), final_ind_classes_size);
//
//   return wrap(class_counts_list);
// }
