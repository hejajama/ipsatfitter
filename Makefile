
include filelist.m

CXX = g++-7
CC = gcc-7
FC = gfortran

MINUITINC = ../inc
MINUITLIBDIR = ../lib/

CXXFLAGS = `gsl-config --cflags` -O3 -I$(MINUITINC) 
FFLAGS = -O2
LDFLAGS = `gsl-config --libs` -lgfortran

all: ipsatfit

ipsatfit: $(OBJECTS) $(COBJECTS) $(FOBJECTS) $(SOURCES) $(CSOURCES)
	$(CXX)  $(CXXFLAGS) $(LDFLAGS) $(MINUITLIBDIR)/libMinuit2.a  $(OBJECTS) $(COBJECTS) $(FOBJECTS) -o ipsatfit

dipole: src/dipoleamplitude.o src/ugd_from_ipsat.o src/dglap_cpp/AlphaStrong.o src/dglap_cpp/EvolutionLO.o src/dglap_cpp/EvolutionLO_nocoupling.o src/ugd_from_ipsat.o 
	$(CXX)  $(CXXFLAGS) $(LDFLAGS) -o dipoleamplitude src/ugd_from_ipsat.o src/dipoleamplitude.o src/dglap_cpp/AlphaStrong.o src/dglap_cpp/EvolutionLO.o  src/dglap_cpp/EvolutionLO_nocoupling.o

.cpp.o: 
	$(CXX) $(CXXFLAGS) $< -c -o $@

.c.o:
	$(CC) $(CXXFLAGS) $< -c -o $@

.f.o:
	$(FC) $(CXXFLAGS) $(FFLAGS) -c $< -c -o $@


clean:
	rm -f $(OBJECTS)
	rm -f $(COBJECTS)
	rm -f $(FOBJECTS)
	rm -f ipsatfit
