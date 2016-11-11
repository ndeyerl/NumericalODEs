# Makefile for LMM example
#
# Daniel R. Reynolds
# Math 6321 @ SMU
# Fall 2016

CXX      = g++
CXXFLAGS = -O --std=c++11
#CXXFLAGS = -O0 -g --std=c++11

HEADERS = lmm.hpp matrix.hpp trapezoidal.hpp newton.hpp resid.hpp rhs.hpp

# executable targets
all : prob4.exe

prob4.exe : prob4.o trapezoidal.o lmm.o newton.o matrix.o
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o : %.cpp $(HEADERS)
	$(CXX) -c $(CXXFLAGS) $< -o $@

# utilities
clean :
	\rm -rf *.txt *.o *.exe *~
