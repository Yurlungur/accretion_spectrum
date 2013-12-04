# Makefile for the accretion_spectrum package.
# Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
# Time-stamp: <2013-12-04 16:22:28 (jonah)>

# The default compiler is g++
CXX = g++

# Flags for the compiler. Ask for warnings. Enable the debugger.
CXXFLAGS = -Wall -g

default: bisection_test

bisection_test: bisection_test_driver.bin
bisection_test_driver.bin: bisection_test_driver.o rkf45.o bisection.o
	$(CXX) $(CXXFLAGS) -o $@ $^
bisection_test_driver.o: bisection.o rkf45.hpp bisection.hpp
bisection.o: bisection.hpp rkf45.hpp

rkf45.o: rkf45.hpp

clean:
	$(RM) rkf45.o bisection_test_driver.o bisection_test_driver.bin