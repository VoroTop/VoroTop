# research
# VoroTop makefile
#
# Author : Emanuel A. Lazar (UPenn) 
# Email  : mLazar@seas.upenn.edu
# Date   : December 5, 2014

# C++ compiler
CXX=g++ -g

# Flags for the C++ compiler
LFLAGS=-L./voro++-0.5/src -lvoro++ 
E_INC=-I./voro++-0.5/src

vorotop : subsystem vorotop.o wvectors.o variables.o filters.o import.o functions.o output.o mtrand.o
	$(CXX) import.o vorotop.o wvectors.o variables.o filters.o functions.o output.o mtrand.o $(LFLAGS) -o VoroTop

subsystem:
	$(MAKE) -C voro++-0.5

vorotop.o : vorotop.cc
	$(CXX) $(E_INC) -c vorotop.cc

wvectors.o : wvectors.cc
	$(CXX) $(E_INC) -c wvectors.cc

variables.o : variables.cc
	$(CXX) $(E_INC) -c variables.cc

filters.o : filters.cc
	$(CXX) $(E_INC) -c filters.cc

import.o : import.cc
	$(CXX) $(E_INC) -c import.cc

functions.o : functions.cc
	$(CXX) $(E_INC) -c functions.cc

output.o : output.cc
	$(CXX) $(E_INC) -c output.cc

mtrand.o : mtrand.cpp
	$(CXX) $(E_INC) -c mtrand.cpp

zip :
	zip -r package.zip *.cc *.hh *.cpp *.h LICENSE README Makefile voro++-0.5 

