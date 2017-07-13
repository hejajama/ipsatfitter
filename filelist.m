SOURCES = src/main.cpp src/wave_function.cpp  src/virtual_photon.cpp src/ipsat.cpp src/dis.cpp src/data.cpp src/interpolation2d.cpp src/dglap_sartre/Dglap.cpp src/dglap_sartre/DglapEvolution.cpp src/dglap_sartre/Laguerre.cpp
FSOURCES = src/LO_evolution_routine.f src/alphaS.f
OBJECTS=$(SOURCES:.cpp=.o)
FOBJECTS = $(FSOURCES:.f=.o)
