# Makefile for homework 3 problems 5 and 2.
#
# Nicole Deyerl
# Math 6321 @ SMU
# Fall 2016

CXX      = g++ --std=c++11
CXXFLAGS = -O

# executable targets
all : prob2.exe prob5.exe

prob2.exe : prob2.cpp adapt_euler.cpp matrix.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

prob5.exe : prob5.o trapezoidal.o newton.o matrix.o
	$(CXX) $(CXXFLAGS) -o $@ $^

prob5.o : prob5.cpp
	$(CXX) $(CXXFLAGS) -c $^

trapezoidal.o : trapezoidal.cpp
	$(CXX) $(CXXFLAGS) -c $^

newton.o : newton.cpp
	$(CXX) $(CXXFLAGS) -c $^

matrix.o : matrix.cpp
	$(CXX) $(CXXFLAGS) -c $^

# utilities
clean :
	\rm -rf *.txt *.o *.exe *~
