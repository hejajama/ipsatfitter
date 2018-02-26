SOURCES = src/main.cpp src/wave_function.cpp  src/virtual_photon.cpp src/ipsat.cpp src/dis.cpp src/data.cpp src/dglap_cpp/AlphaStrong.cpp src/dglap_cpp/EvolutionLO.cpp src/interpolation.cpp src/woodsaxon.cpp src/dglap_cpp/EvolutionLO_nocoupling.cpp
CSOURCES =   #src/dglap_sartre/dglap.c src/dglap_sartre/laguerre.c
FSOURCES = src/LO_evolution_routine.f src/alphaS.f
OBJECTS=$(SOURCES:.cpp=.o)
COBJECTS=$(CSOURCES:.c=.o)
FOBJECTS = $(FSOURCES:.f=.o)
# src/dglap_sartre/AlphaStrong.cpp
#src/dglap_sartre/DglapEvolution.cpp
