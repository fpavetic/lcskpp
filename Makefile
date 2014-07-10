all:
	g++ -o test_lcskpp test_lcskpp.cpp lcskpp.cpp -O2 -std=c++11

clean:
	rm -f test_lcskpp

