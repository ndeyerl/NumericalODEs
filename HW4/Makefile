# Makefile for homework 2
#
# Daniel R. Reynolds
# Math 6321 @ SMU
# Fall 2016

# compiler & flags
CXX = g++
CXXFLAGS = -O -std=c++11
#CXXFLAGS = -O0 -g -std=c++11

# executable targets
all: prob2.exe prob4.exe 

prob2.exe : prob2.cpp erk4.cpp ab3.cpp matrix.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

prob4.exe : prob4.cpp adapt_rkf.cpp adapt_euler.cpp matrix.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

# utilities
clean :
	\rm -rf *.txt *.png *.exe *~
