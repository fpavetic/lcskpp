/**
 * @author: Filip Pavetic (fpavetic@gmail.com)
 */

#include <algorithm>
#include <map>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#include <cassert>
#include <cmath>
#include <cstdlib>

#include "fenwick.h"
#include "lcskpp.h"
using namespace std;

// This function determines the total number of 
// distinct characters in input strings a and b.
// Outputs are: aid[character] = unique_character_id
//              alphabet_size = total number of distinct chars
static void prepare_alphabet(
    const string& a, const string& b,
    vector<char>& aid, int& alphabet_size) {
  aid = vector<char>(256, -1);
  alphabet_size = 0;

  for (size_t i = 0; i < a.size(); ++i) {
    if (aid[a[i]] == -1) {
      aid[a[i]] = alphabet_size++;
    }
  }
  for (size_t i = 0; i < b.size(); ++i) {
    if (aid[b[i]] == -1) {
      aid[b[i]] = alphabet_size++;
    }
  }
}

// matches is filled with sorted pairs (i,j) meaning that for 
// every such pair a[i...i+k-1] == b[j...j+k-1].
//
// In this implementation it is assumed that 
// power(alphabet_size, k) < power(2, 64). In case this is
// not true in the particular application, this function
// would have to be reimplemented (for example by using
// non-perfect hashing or suffix arrays) and everything else
// should work again.
static void get_matches(
    const string& a, const string& b,
    const int k, vector<pair<int, int> >* matches) {
  assert(matches != NULL);
  matches->clear();

  vector<char> aid;
  int alphabet_size;
  prepare_alphabet(a, b, aid, alphabet_size);

  // We assume: alphabet_size ** k < 2 ** 64,
  // in the case this does not hold, the entire
  // get_matches function probably has to be
  // reimplemented using non-perfect hashing or
  // suffix arrays.
  if (k * log(alphabet_size) >= 64 * log(2)) {
    fprintf(stderr, "We assume that alphabet_size ** k <\
2 ** 64.\nPlease see lcskpp.cpp for more information.");
    exit(1);
  }

  typedef unordered_multimap<uint64_t, int> MatchIndexType;
  unique_ptr<MatchIndexType> match_index =
    unique_ptr<MatchIndexType>(new MatchIndexType());

  uint64_t hash_mod = 1;
  for (int i = 0; i < k; ++i) hash_mod *= alphabet_size;

  if (alphabet_size == 4) {
    assert(hash_mod == (1LL<<(2*k)));
  }

  uint64_t rolling_hash = 0;
  for (int i = 0; i < a.size(); ++i) {
    rolling_hash = rolling_hash * alphabet_size + aid[a[i]];
    rolling_hash %= hash_mod;

    if (i+1 >= k) {
      match_index->insert(
        MatchIndexType::value_type(rolling_hash, i-k+1));
    }
  }

  rolling_hash = 0;
  for (int i = 0; i < b.size(); ++i) {
    rolling_hash = rolling_hash * alphabet_size + aid[b[i]];
    rolling_hash %= hash_mod;

    if (i+1 >= k) {
      auto positions_in_a = match_index->equal_range(rolling_hash);
      for (auto it = positions_in_a.first; 
           it != positions_in_a.second; ++it) {
        matches->push_back(make_pair(it->second, i-k+1));
      }
    }
  }

  sort(matches->begin(), matches->end());
}

// This function is used to reconstruct the indices
// of the calculated LCSk++. lcskpp_recon is going to
// contain pairs of indices of matching characters.
static void fill_lcskpp_reconstruction(
    const vector<pair<int, int> >& matches,
    const int k,
    const vector<int>& prev_idx,
    const int last_idx,
    const int lcskpp_len,
    vector<pair<int, int> >* lcskpp_recon) {
  assert(lcskpp_recon != NULL);
  lcskpp_recon->clear();
  lcskpp_recon->reserve(lcskpp_len);

  for (int i = last_idx; i != -1; i = prev_idx[i]) {
    int r = matches[i].first+k-1;
    int c = matches[i].second+k-1;

    if (prev_idx[i] == -1 || 
        (matches[prev_idx[i]].first + k <= matches[i].first &&
         matches[prev_idx[i]].second + k <= matches[i].second)) {
      // Taking the entire match ...
      for (int j = 0; j < k; ++j, --r, --c) {
        lcskpp_recon->push_back(make_pair(r, c));
      }
    } else { 
      // ... otherwise it is a continuation. Only the
      // difference is taken.
      int curr = i;
      int prev = prev_idx[i];

      int curr_secondary_diag = 
        (matches[curr].first + matches[curr].second) / 2;
      int prev_secondary_diag = 
        (matches[prev].first + matches[prev].second) / 2;
      assert(prev_secondary_diag < curr_secondary_diag);

      for (int j = prev_secondary_diag; 
           j < curr_secondary_diag; ++j, --r, --c) {
        lcskpp_recon->push_back(make_pair(r, c));
      }
    }
  }

  reverse(lcskpp_recon->begin(), lcskpp_recon->end());
}

// Assumes matches is sorted by standard pair ordering.
static void lcskpp_sparse_fast(
    const vector<pair<int, int> >& matches,
    const int k, int* lcskpp_length,
    vector<pair<int, int> >* lcskpp_reconstruction) {
  if (matches.empty()) {
    *lcskpp_length = 0;
    if (lcskpp_reconstruction) lcskpp_reconstruction->clear();
    return;
  }

  vector<tuple<int, int, int> > events;
  events.reserve(2*matches.size());

  int n = 0;
  for (auto it = matches.begin(); it != matches.end(); ++it) {
    int idx = it - matches.begin();
    events.push_back(make_tuple(it->first, it->second, 
				idx+matches.size())); // begin
    events.push_back(make_tuple(it->first+k, it->second+k, idx)); // end
    
    n = max(n, it->first+k);
    n = max(n, it->second+k);
  }
  sort(events.begin(), events.end());

  // Indexed by column, first:dp value, second:index in matches.
  FenwickMax<pair<int, int> > dp_col_max(n);
  vector<int> dp(matches.size());
  vector<int> recon(matches.size());
  vector<int> continues(matches.size(), -1);
  if (k > 1) {
    for (auto curr = matches.begin(); 
         curr != matches.end(); ++curr) {
      auto G = make_pair(curr->first-1, curr->second-1);
      auto prev = lower_bound(matches.begin(), matches.end(), G);
      if (*prev == G) {
	continues[curr-matches.begin()] = prev-matches.begin();
      }
    }
  }

  int best_idx = 0;
  *lcskpp_length = 0;
    
  for (auto event = events.begin(); 
       event != events.end(); ++event) {
    int idx = get<2>(*event) % matches.size();
    bool is_beginning = (get<2>(*event) >= matches.size());
    int i = get<0>(*event);
    int j = get<1>(*event);
    int primary_diagonal = n-1+i-j;

    if (is_beginning) { // begin
      pair<int, int> prev_dp = dp_col_max.get(j);
      dp[idx] = k;
      recon[idx] = -1;

      if (prev_dp.first > 0) {
	dp[idx] = prev_dp.first + k;
	recon[idx] = prev_dp.second;
      }
    } else {
      if (continues[idx] != -1) {
	if (dp[continues[idx]] + 1 > dp[idx]) {
	  dp[idx] = dp[continues[idx]] + 1;
	  recon[idx] = continues[idx];
	}
      }

      dp_col_max.update(j, make_pair(dp[idx], idx));

      if (dp[idx] > *lcskpp_length) {
	*lcskpp_length = dp[idx];
	best_idx = idx;
      }
    }
  }
  
  if (lcskpp_reconstruction) {
    fill_lcskpp_reconstruction(matches, k, recon, best_idx, 
                               *lcskpp_length,
                               lcskpp_reconstruction);
  }    
}

// Assume matches is sorted by standard pair ordering.
static void lcskpp_sparse_slow(
    const vector<pair<int, int> >& matches,
    const int k, int* lcskpp_length,
    vector<pair<int, int> >* lcskpp_reconstruction) {
  assert(is_sorted(matches.begin(), matches.end()));

  if (matches.empty()) {
    *lcskpp_length = 0;
    if (lcskpp_reconstruction) lcskpp_reconstruction->clear();
  } else {
    int n = matches.size();
    vector<int> dp(n);
    vector<int> recon(n);
    int best_idx = 0;
    *lcskpp_length = 0;
    
    for (int i = 0; i < n; ++i) {
      dp[i] = k;
      recon[i] = -1;

      int end_row_i = matches[i].first+k-1;
      int end_col_i = matches[i].second+k-1;

      int primary_diagonal_i = end_row_i - end_col_i;
      int secondary_diagonal_i = (end_row_i + end_col_i)/2;

      for (int j = i-1; j >= 0; --j) {
        if (matches[j].first + k <= matches[i].first &&
            matches[j].second + k <= matches[i].second) {
          // 1) Taking the entire match interval and continuing
          // another match which 'ended'.
          if (dp[j] + k > dp[i]) {
            dp[i] = dp[j] + k;
            recon[i] = j;
          }
        } else {
          // 2) Continuation on the same primary diagonal.
          int end_row_j = matches[j].first+k-1;
          int end_col_j = matches[j].second+k-1;
          int primary_diagonal_j = end_row_j - end_col_j;
          int secondary_diagonal_j = (end_row_j + end_col_j)/2;
	  
          if (primary_diagonal_i == primary_diagonal_j &&
              secondary_diagonal_i > secondary_diagonal_j &&
              secondary_diagonal_i - secondary_diagonal_j < k) {
            if (dp[j] + secondary_diagonal_i - secondary_diagonal_j 
                > dp[i]) {
              dp[i] = 
                dp[j] + secondary_diagonal_i - secondary_diagonal_j;
              recon[i] = j;
            }
          }
        }
      }

      if (dp[i] > *lcskpp_length) {
        best_idx = i;
        *lcskpp_length = dp[i];
      }
    }

    if (lcskpp_reconstruction) {
      fill_lcskpp_reconstruction(matches, k, recon, best_idx, 
                                 *lcskpp_length, 
                                 lcskpp_reconstruction);
    }
  }
}

void lcskpp_sparse_slow(
    const string& a, const string& b,
    const int k, int* lcskpp_length,
    vector<pair<int, int> >* lcskpp_reconstruction) {
  vector<pair<int, int> > matches;
  get_matches(a, b, k, &matches);
  lcskpp_sparse_slow(matches, k, lcskpp_length, 
                     lcskpp_reconstruction);  
}

void lcskpp_sparse_fast(
    const string& a, const string& b,
    const int k, int* lcskpp_length,
    vector<pair<int, int> >* lcskpp_reconstruction) {
  vector<pair<int, int> > matches;
  get_matches(a, b, k, &matches);
  lcskpp_sparse_fast(matches, k, lcskpp_length, 
                     lcskpp_reconstruction);  
}

// Assume matches is sorted by standard pair ordering.
void lcskpp(const vector<pair<int, int> >& matches,
          const int k, int* lcskpp_length,
          vector<pair<int, int> >* lcskpp_reconstruction) {
  if (matches.size() < 700) {
    lcskpp_sparse_slow(matches, k, lcskpp_length, 
                       lcskpp_reconstruction);
  } else {
    lcskpp_sparse_fast(matches, k, lcskpp_length, 
                       lcskpp_reconstruction);
  }
}

void lcskpp(const string& a, const string& b,
          const int k, int* lcskpp_length,
          vector<pair<int, int> >* lcskpp_reconstruction) {
  vector<pair<int, int> > matches;
  get_matches(a, b, k, &matches);
  lcskpp(matches, k, lcskpp_length, lcskpp_reconstruction);  
}

bool valid_lcskpp(const string& a, const string& b,
                  const int k, const int lcskpp_len,
                  const vector<pair<int, int> >& lcskpp_recon) {
  // 1) Ensure correct length.
  if (lcskpp_len != lcskpp_recon.size()) {
    return false;
  }

  // 2) Ensure chars corresponding to the indices match.
  for (auto match: lcskpp_recon) {
    int i = match.first;
    int j = match.second;
    
    if (i < 0 || i >= a.size()) { return false; }
    if (j < 0 || j >= b.size()) { return false; }
    if (a[i] != b[j]) { return false; }
  }

  // 3) Ensure runs of indices have at least length of k.
  int run_a = 1;
  int run_b = 1;
  for (size_t i = 1; i < lcskpp_recon.size(); ++i) {
    if (lcskpp_recon[i-1].first >= lcskpp_recon[i].first) { return false; }
    if (lcskpp_recon[i-1].second >= lcskpp_recon[i].second) { return false; }

    if (lcskpp_recon[i-1].first+1 == lcskpp_recon[i].first) { ++run_a; }
    if (lcskpp_recon[i-1].second+1 == lcskpp_recon[i].second) { ++run_b; }

    if (i+1 == lcskpp_recon.size() || 
        lcskpp_recon[i-1].first+1 != lcskpp_recon[i].first) {
      if (run_a < k) { return false; }
      run_a = 1;
    }

    if (i+1 == lcskpp_recon.size() ||
        lcskpp_recon[i-1].second+1 != lcskpp_recon[i].second) {
      if (run_b < k) { return false; }
      run_b = 1;
    }
  }

  return true;
}

static int min3(int a, int b, int c) { return min(min(a,b),c); }

void lcskpp_slow(const string& a, const string& b, const int K,
                 int* lcskpp_length) {
  vector<vector<int> > dp(a.size()+1, 
                          vector<int>(b.size()+1));

  for (int i = 0; i <= a.size(); ++i) dp[i][0] = 0;
  for (int j = 0; j <= b.size(); ++j) dp[0][j] = 0;

  for (int i = 1; i <= a.size(); ++i) {
    for (int j = 1; j <= b.size(); ++j) {
      dp[i][j] = max(dp[i-1][j], dp[i][j-1]);

      // 2*K is good enough limit because everything bigger is 
      // covered by some precedent matching interval.
      for (int k = 1; k <= min3(i, j, 2*K); ++k) {
        char aa = a[i-k];
        char bb = b[j-k];
        if (aa != bb) {
          break;
        }

        if (k >= K) {
          dp[i][j] = max(dp[i][j], dp[i-k][j-k]+k);
        }
      }
    }
  }

  *lcskpp_length = dp[a.size()][b.size()];
}
