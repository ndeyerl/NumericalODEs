# Makefile for forward Euler examples
#
# Daniel R. Reynolds
# SMU Mathematics
# Fall 2016

# compiler & flags
CXX = g++
CXXFLAGS = -O -std=c++11
#CXXFLAGS = -O0 -g -std=c++11


# makefile targets
all : prob1.exe prob3.exe driver_fwd_euler.exe driver_fwd_euler_system.exe

prob1.exe : prob1.cpp fwd_euler.cpp matrix.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^
	
prob3.exe : prob3.cpp adapt_euler.cpp matrix.cpp 
	$(CXX) $(CXXFLAGS) -o $@ $^

driver_fwd_euler.exe : driver_fwd_euler.cpp fwd_euler.cpp matrix.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

driver_fwd_euler_system.exe : driver_fwd_euler_system.cpp fwd_euler.cpp matrix.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^



# utilities
clean :
	\rm -rf *.txt *.exe *~
