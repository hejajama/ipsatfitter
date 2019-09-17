
include filelist.m

CXX = g++-8
CC = gcc-8
FC = gfortran

MINUITINC = ../minuit2//include/
MINUITLIBDIR = ../minuit2/lib/

CXXFLAGS = `gsl-config --cflags` -O3 -I$(MINUITINC) 
FFLAGS = -O2
LDFLAGS = `gsl-config --libs` -lgfortran

all: ipsatfit

ipsatfit: src/main.o $(OBJECTS) $(COBJECTS) $(FOBJECTS) $(SOURCES) $(CSOURCES)
	$(CXX)  $(CXXFLAGS) $(LDFLAGS) $(MINUITLIBDIR)/libMinuit2.a  src/main.o $(OBJECTS) $(COBJECTS) $(FOBJECTS) -o ipsatfit

dipole: src/dipoleamplitude.o src/ugd_from_ipsat.o src/dglap_cpp/AlphaStrong.o  src/dglap_cpp/EvolutionLO_nocoupling.o src/ugd_from_ipsat.o 
	$(CXX)  $(CXXFLAGS) $(LDFLAGS) -o dipoleamplitude src/ugd_from_ipsat.o src/dipoleamplitude.o src/dglap_cpp/AlphaStrong.o  src/dglap_cpp/EvolutionLO_nocoupling.o

tools: $(OBJECTS) tools/create_light_f2.o 
	$(CXX)  $(CXXFLAGS) $(LDFLAGS) -o tools/create_light_f2 $(MINUITLIBDIR)/libMinuit2.a tools/create_light_f2.o $(OBJECTS) $(COBJECTS) $(FOBJECTS)

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
