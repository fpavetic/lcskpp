/**
 * Utilities to generate random strings as in the
 * similarity model described in the paper.
 * 
 * @author Filip Pavetic (fpavetic@gmail.com)
 */

#include <cassert>
#include <cstdlib>
#include <string>

const std::string kNuc = "ACTG";

// Selecting a character uniformily.
char get_random_base() {
  return kNuc[rand()%4];
}

// This function returns a string which is a copy
// the input string, but character on every position
// is mutated to another one from the alphabet
std::string generate_similar(const std::string& a, 
                             const double& p_err) {
  std::string b = a;
  for (int i = 0; i < b.size(); ++i) {
    if (1.0*rand()/RAND_MAX <= p_err) {
      b[i] = get_random_base();
    }
  }
  return b;
}

// This generates a random string of length <len>.
std::string generate_string(const int& len) {
  std::string ret;
  for (int i = 0; i < len; ++i) {
    ret += get_random_base();
  }
  assert(ret.size() == len);
  return ret;
}
