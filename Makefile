
include filelist.m

CXX = /opt/local/bin/g++-mp-5

MINUITINC = ../inc
MINUITLIBDIR = ../lib/

CXXFLAGS = `/opt/local/bin/gsl-config --cflags` -O2 -I$(MINUITINC) -fopenmp
LDFLAGS = `/opt/local/bin/gsl-config --libs` -lgfortran
all: ipsatfit

ipsatfit: $(OBJECTS)
	$(CXX)  $(CXXFLAGS) $(LDFLAGS) $(MINUITLIBDIR)/libMinuit2.a $(OBJECTS) libColorDipole/libraries/libColorDipole.a -o ipsatfit

.cpp.o: 
	$(CXX) $(CXXFLAGS) $< -c -o $@

clean:
	rm -f $(OBJECTS)
	rm -f ipsatfit
