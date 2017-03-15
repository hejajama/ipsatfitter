
include filelist.m

CXX = /opt/local/bin/g++-mp-5
FC = /opt/local/bin/gfortran-mp-5

MINUITINC = ../inc
MINUITLIBDIR = ../lib/

CXXFLAGS = `/opt/local/bin/gsl-config --cflags` -O2 -I$(MINUITINC) -fopenmp
FFLAGS = -O2
LDFLAGS = `/opt/local/bin/gsl-config --libs` -lgfortran

all: ipsatfit

ipsatfit: $(OBJECTS) $(FOBJECTS)
	$(CXX)  $(CXXFLAGS) $(LDFLAGS) $(MINUITLIBDIR)/libMinuit2.a $(OBJECTS) $(FOBJECTS)  -o ipsatfit

.cpp.o: 
	$(CXX) $(CXXFLAGS) $< -c -o $@

.f.o:
	$(FC) $(FFLAGS) -c $< -c -o $@


clean:
	rm -f $(OBJECTS)
	rm -f $(FOBJECTS)
	rm -f ipsatfit
