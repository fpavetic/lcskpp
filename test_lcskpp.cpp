/**
 * This is a tool used to test the implementation of the
 * algorithms for calculating the LCSk++.
 *
 * @author: Filip Pavetic (fpavetic@gmail.com)
 */

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "lcskpp.h"
#include "random_strings.h"
using namespace std;

// Length of the strings.
const int kStringLen = 100; 

// Number of performed simulations.
const int kSimulationRuns = 1000; 

// Default value of the k parameter.
const int kK = 3;

// If kPerr is set to -1 then rand-to-rand strings are aligned, otherwise rand-to-modified-copy simulations are performed.
// const double kPerr = -1.0;
const double kPerr = 0.05;

int test_lcskpp(const string& a, const string& b, const int K) {
  int lcskpp_length = 0;
  lcskpp_slow(a, b, K, &lcskpp_length);

  int lcskpp_sparse_slow_len = 0;
  vector<pair<int, int> > lcskpp_sparse_slow_recon;
  lcskpp_sparse_slow(a, b, K, &lcskpp_sparse_slow_len,
                     &lcskpp_sparse_slow_recon);

  int lcskpp_sparse_fast_len = 0;
  vector<pair<int, int> > lcskpp_sparse_fast_recon;
  lcskpp_sparse_fast(a, b, K, &lcskpp_sparse_fast_len,
                     &lcskpp_sparse_fast_recon);
  
  // printf("%d %d %d\n", 
  //        lcskpp_length, 
  //        lcskpp_sparse_slow_len,
  //        lcskpp_sparse_fast_len);
  // printf("%s\n%s\n", a.c_str(), b.c_str());
  assert(lcskpp_length == lcskpp_sparse_slow_len);
  assert(lcskpp_length == lcskpp_sparse_fast_len);
  assert(valid_lcskpp(a, b, K, lcskpp_sparse_slow_len, 
		      lcskpp_sparse_slow_recon));
  assert(valid_lcskpp(a, b, K, lcskpp_sparse_fast_len,
    		      lcskpp_sparse_fast_recon));
  
  return lcskpp_length;
}

int run_one_simulation() {
  pair<string, string> ab;
  ab.first = generate_string(kStringLen);
  ab.second = kPerr < 0 ? 
                      generate_string(kStringLen) : 
                      generate_similar(ab.first, kPerr);
  const string& a = ab.first;
  const string& b = ab.second;
  return test_lcskpp(a, b, kK);
}


void calculate_distribution(map<int, double>& distr) {
  distr.clear();
  for (int i = 0; i < kSimulationRuns; ++i) {
    distr[run_one_simulation()] += 1.0/kSimulationRuns;
  }
}

int main(int argc, char* argv[]) {
  srand(1603);

  map<int, double> distr;
  calculate_distribution(distr);

  double sum_prob = 0;
  double e_lcs = 0;

  for (int i = 0; i <= kStringLen; ++i) {
    double p = distr[i];
    sum_prob += p;
    e_lcs += p*i;
  }

  assert(0.99999 <= sum_prob <= 1.00001);
  printf("Expected k-LCS=%0.3lf\n", e_lcs);
  return 0;
}
