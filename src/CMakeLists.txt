configure_file("${CMAKE_CURRENT_SOURCE_DIR}/gitsha1.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp" @ONLY)

add_library(
	ipsatdglap
	STATIC
	dglap_cpp/AlphaStrong.cpp
	dglap_cpp/EvolutionLO_nocoupling.cpp
	)

add_executable(
	ipsatfit	
	main.cpp
	ipsat.cpp
	dis.cpp
	data.cpp
	woodsaxon.cpp
	interpolation.cpp
	virtual_photon.cpp
	wave_function.cpp
	${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp
	LO_evolution_routine.f 
	alphaS.f
	)

include_directories(.)

find_package(GSL REQUIRED)    #

target_link_libraries(
	ipsatfit
	PRIVATE
	Minuit2.a
	ipsatdglap
	GSL::gsl
	GSL::gslcblas
)


