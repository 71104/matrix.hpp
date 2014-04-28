CXXFLAGS+=-Wall --std=c++11 -Ilib

default : bin/test

bin :
	mkdir -p bin

bin/test : bin src/test.cpp lib/matrix.hpp
	$(CXX) $(CXXFLAGS) -o bin/test src/test.cpp

clean :
	rm -rf ./bin
