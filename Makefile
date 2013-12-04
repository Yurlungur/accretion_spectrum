# Makefile for the accretion_spectrum package.
# Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
# Time-stamp: <2013-12-04 18:16:23 (jonah)>

# The default compiler is g++
CXX = g++

# Flags for the compiler. Ask for warnings. Enable the debugger.
CXXFLAGS = -Wall -g

default: accretion_spectrum_test
all: bisection_test,accretion_spectrum_test

accretion_spectrum_test: accretion_spectrum_test_driver.bin
accretion_spectrum_test_driver.bin: accretion_spectrum_test_driver.o accretion_spectrum.o bisection.o rkf45.o
	$(CXX) $(CXXFLAGS) -o $@ $^

accretion_spectrum_test_driver.o: accretion_spectrum.hpp bisection.hpp rkf45.hpp
accretion_spectrum.o: accretion_spectrum.hpp bisection.hpp rkf45.hpp

bisection_test: bisection_test_driver.bin
bisection_test_driver.bin: bisection_test_driver.o rkf45.o bisection.o
	$(CXX) $(CXXFLAGS) -o $@ $^
bisection_test_driver.o: bisection.o rkf45.hpp bisection.hpp
bisection.o: bisection.hpp rkf45.hpp

rkf45.o: rkf45.hpp

.PHONY: default, all, bisection_test, accretion_spectrum_test

clean:
	$(RM) rkf45.o bisection_test_driver.o bisection_test_driver.bin