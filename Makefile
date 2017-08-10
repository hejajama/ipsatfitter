
include filelist.m

CXX = /opt/local/bin/g++-mp-5
CC = /opt/local/bin/gcc-mp-5
FC = /opt/local/bin/gfortran-mp-5

MINUITINC = ../inc
MINUITLIBDIR = ../lib/

CXXFLAGS = `/opt/local/bin/gsl-config --cflags` -g -I$(MINUITINC)
FFLAGS = -O2
LDFLAGS = `/opt/local/bin/gsl-config --libs` -lgfortran

all: ipsatfit

ipsatfit: $(SOURCES) $(CSOURCES)
	$(CXX)  $(CXXFLAGS) $(LDFLAGS) $(MINUITLIBDIR)/libMinuit2.a  $(OBJECTS) $(COBJECTS) $(FOBJECTS)  -o ipsatfit

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
