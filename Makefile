CXXFLAGS+=-Wall --std=c++11 -I.

test : test.cpp matrix.hpp
	$(CXX) $(CXXFLAGS) -o test test.cpp

clean :
	rm -rf ./test
