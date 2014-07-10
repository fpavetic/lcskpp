LCSk++: Practical similarity metric for long strings
======

This is an implementation of the LCSk++ metric for long strings described in [1].

LCSk++ of two strings a and b calculates the longest common subsequence of two strings with
the restriction that the consecutive runs of indices in both strings have length at least k,
which is a parameter of the algorithm. For example: longest common subsequence of the strings
ABCDAB and ABCADB is of length 5 (ABCDB), LCSk++ of these two strings with k=3 is 3 (ABC).
This restriction loses some matches, but allows for a faster computation.

Implementation
======
lcskpp.h/lcskpp.cpp: Implementation of several algorithms for computing LCSk++. The approach described 
                     in Section 3.2 of [1] can be found in the lcskpp_sparse_fast function.
fenwick.h: Implementation of the Fenwick tree data structure used by the lcskpp_sparse_fast.                     
test_klcs.cpp: A unit test for the algorithm.
random_strings.h: Functions for generating random strings as described in Section 4.1 of [1].
                     
Dependencies
============
For compiling this library, it is necessary to have C++11 compatible compiler.

References
==========
[1] Filip Pavetic, Goran Zuzic, Mile Sikic: "LCSk++: Practical similarity metric for long strings" 
    (http://arxiv.org/pdf/1407.2407v1.pdf)
    This paper has been created as a continuation of the Master Thesis of the first author,
    made on the Faculty of Electrical Engineering and Computing at University of Zagreb, Croatia
