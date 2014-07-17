### LCSk++: Practical similarity metric for long strings

This is an implementation of the LCSk++ metric for long strings described in [1].

LCSk++ of two strings a and b calculates the longest common subsequence of two strings with
the restriction that the consecutive runs of indices in both strings have length at least k,
which is a parameter of the algorithm. For example: longest common subsequence of the strings
ABCDAB and ABCADB is of length 5 (ABCDB), while LCSk++ of these two strings with k=3 is 3 (ABC).
This restriction loses some matches, but allows for a faster computation.

### Implementation
* __lcskpp.h/lcskpp.cpp__  
   >> Implementation of several algorithms for computing LCSk++. The approach described in Section 3.2 of [1] can be found in the lcskpp_sparse_fast function.
* __fenwick.h__  
   >> Implementation of the Fenwick tree data structure used by the lcskpp_sparse_fast.
* __test_lcskpp.cpp__  
   >> A unit test for the algorithm.
* __random_strings.h__  
   >> Functions for generating random strings as described in Section 4.1 of [1].

### Dependencies
For compiling this library, it is necessary to have C++11 compatible compiler.

### References
[1] Filip Pavetic, Goran Zuzic, Mile Sikic: _LCSk++: Practical similarity metric for long strings_, http://arxiv.org/abs/1407.2407  
[2] Gary Benson, Avivit Levy, Riva Shalom: _Longest Common Subsequence in k-length substrings_, http://arxiv.org/abs/1402.2097  
[3] Sebastian Deorowicz, Szymon Grabowski: _Efficient algorithms for the longest common subsequence in k-length substrings_, http://arxiv.org/abs/1311.4552

Note: [1] has been created as a continuation of the first authors Master Thesis, written on the Faculty of Electrical Engineering and Computing, University of Zagreb, Croatia.
