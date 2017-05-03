# research
# VoroTop makefile
#
# Author : Emanuel A. Lazar (UPenn) 
# Email  : mLazar@seas.upenn.edu
# Date   : December 5, 2014

# C++ compiler
CXX      = g++ 

# Flags for the C++ compiler
LFLAGS   = -L./voro++-0.5/src -lvoro++ 
LIBS     = -I./voro++-0.5/src
CXXFLAGS = -Wall -O3 -c -std=c++11 


vorotop : subsystem vorotop.o wvectors.o variables.o filters.o import.o functions.o output.o 
	$(CXX) import.o vorotop.o wvectors.o variables.o filters.o functions.o output.o $(LFLAGS) -o VoroTop

subsystem:
	$(MAKE) -C voro++-0.5

vorotop.o : vorotop.cc import.hh filters.hh vorotop.hh functions.hh
	$(CXX) $(CXXFLAGS) $(LIBS) vorotop.cc

wvectors.o : wvectors.cc filters.hh vorotop.hh
	$(CXX) $(CXXFLAGS) $(LIBS) wvectors.cc

variables.o : variables.cc
	$(CXX) $(CXXFLAGS) $(LIBS) variables.cc

filters.o : filters.cc filters.hh vorotop.hh
	$(CXX) $(CXXFLAGS) $(LIBS) filters.cc

import.o : import.cc import.hh filters.hh vorotop.hh
	$(CXX) $(CXXFLAGS) $(LIBS) import.cc

functions.o : functions.cc filters.hh vorotop.hh
	$(CXX) $(CXXFLAGS) $(LIBS) functions.cc

output.o : output.cc import.hh filters.hh vorotop.hh
	$(CXX) $(CXXFLAGS) $(LIBS) output.cc

zip :
	zip -r package.zip *.cc *.hh LICENSE README Makefile voro++-0.5 -x voro++-0.5/src/*.o voro++-0.5/src/libvoro++.a voro++-0.5/src/voro++ 

clean :
	rm -f *.o voro++-0.5/src/*.o voro++-0.5/src/libvoro++.a voro++-0.5/src/voro++

