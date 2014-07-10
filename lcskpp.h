/**
 * Implementation of the algorithms for calculating LCSk++ 
 * similarity metric described in the paper:
 * LCSk++: A practical similarity metric for long strings
 * (http://arxiv.org/pdf/1407.2407v1.pdf) 
 *
 * @author: Filip Pavetic (fpavetic@gmail.com)
 */

#ifndef LCSKPP
#define LCSKPP

#include <string>
#include <utility>
#include <vector>

// A calculation of LCSk++ which switches between
// lcskpp_sparse_slow and lcskpp_sparse_fast, 
// depending on the number of match pairs in the
// strings.
// 
// If lcskpp_reconstruction equals NULL, only the
// value of the metric is computed, without reconstructing it.
void lcskpp(
    const std::string& a, const std::string& b, 
    const int k, int* lcskpp_length,
    std::vector<std::pair<int, int> >* lcskpp_reconstruction);

// This is a slower, match pair based, LCSk++ calculation.
// When number of match pairs is low, this can perform better in
// practice.
// 
// If lcskpp_reconstruction equals NULL, only the
// value of the metric is computed, without reconstructing it.
void lcskpp_sparse_slow(
    const std::string& a, const std::string& b, 
    const int k, int* lcskpp_length,
    std::vector<std::pair<int, int> >* lcskpp_reconstruction);

// This is the algorithm described in Section 3.2 of the paper 
// LCSk++: A practical similarity metric for long strings
// (http://arxiv.org/pdf/1407.2407v1.pdf)
// 
// If lcskpp_reconstruction equals NULL, only the
// value of the metric is computed, without reconstructing it.
void lcskpp_sparse_fast(
    const std::string& a, const std::string& b, 
    const int k, int* lcskpp_length,
    std::vector<std::pair<int, int> >* lcskpp_reconstruction);

// This function tests whether LCSk++ has been reconstructed
// successfully.
bool valid_lcskpp(
    const std::string& a, const std::string& b,
    const int k, const int lcskpp_len,
    const std::vector<std::pair<int, int> >& lcskpp_recon);

// Slow calculation of LCSk++, for testing purposes.
void lcskpp_slow(
    const std::string& a, const std::string& b, const int K,
    int* lcskpp_length);

#endif
