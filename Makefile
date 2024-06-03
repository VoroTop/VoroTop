########################################################
####                                                ####
####   ******************************************   ####
####   *                                        *   ####
####   *     VoroTop: Voronoi Cell Topology     *   ####
####   *   Visualization and Analysis Toolkit   *   ####
####   *             (Version 1.0)              *   ####
####   *                                        *   ####
####   *           Emanuel A. Lazar             *   ####
####   *          Bar Ilan University           *   ####
####   *            September 2022              *   ####
####   *                                        *   ####
####   ******************************************   ####
####                                                ####
########################################################

####   File: Makefile

# C++ compiler
CXX      = g++ -std=c++11 -fopenmp

# Flags for the C++ compiler
LFLAGS   = -lvoro++
CXXFLAGS = -Wall -O3 -c


vorotop : vorotop.o vectors.o variables.o filters.o import.o functions.o analysis.o output.o 
	$(CXX) import.o vorotop.o vectors.o variables.o filters.o functions.o analysis.o output.o $(LFLAGS) -o VoroTop

vorotop.o : vorotop.cc filters.hh variables.hh functions.hh
	$(CXX) $(CXXFLAGS) vorotop.cc

vectors.o : vectors.cc filters.hh variables.hh 
	$(CXX) $(CXXFLAGS) vectors.cc

variables.o : variables.cc
	$(CXX) $(CXXFLAGS) variables.cc

filters.o : filters.cc filters.hh variables.hh 
	$(CXX) $(CXXFLAGS) filters.cc

import.o : import.cc filters.hh variables.hh 
	$(CXX) $(CXXFLAGS) import.cc

functions.o : functions.cc filters.hh variables.hh 
	$(CXX) $(CXXFLAGS) functions.cc

analysis.o : analysis.cc filters.hh variables.hh
	$(CXX) $(CXXFLAGS) analysis.cc

output.o : output.cc filters.hh variables.hh 
	$(CXX) $(CXXFLAGS) output.cc

zip :
	zip -r package.zip *.cc *.hh LICENSE README Makefile 

clean :
	rm -f *.o 


