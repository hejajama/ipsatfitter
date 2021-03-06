IPsat fitter to combined HERA data using MINUIT

Reference: H. Mantysaari, P. Zurita, Phys.Rev. D98 (2018) 036002, arXiv:1804.05311

Questions and comments: heikki.mantysaari@jyu.fi

Install:

0. Download MINUIT2 from http://seal.web.cern.ch/seal/snapshot/work-packages/mathlibs/minuit/ -> Download
   I have used Minut2-5.34.14 from
1. Install Minuit2 (CC=compiler CXX=compiler ./configure --prefix=/installation_directory; make; make install)
    Example config:
    CC=/opt/local/bin/gcc-mp-5 CXX=/opt/local/bin/g++-mp-5 ./configure --preix=/Users/heikkimantysaari/code/minuit2/ --disable-openmp
2. Modify Makefile so that C++ and Fortran compiler paths and Minuit paths are correct
3. Compile fitter: make

Parallerization note: if may be a good idea to compile MINUIT with --disable-openmp, and then parallerize chi^2 calculation in DISFitter class. At least if the number of parameters we are fitting is small. Because if MINUIT takes care of parallerization, it calculates chi^2 in parallel with different parametrizations, but with n parameters it ends up running on n threads (experimental fact, if I have misunderstood, then think this more!). If parallerization is done at that level, #pragma omp parallel for in DISFitter can not create any more threads, and thus does not actually parallerize the slow chi^2 calculation!

    When Chi^2 calculation is parallerized, we can calculate in parallel all (x,Q^2) points and
    could run up to ~50 CPUs quite efficiently!


Workflow:
Read datafiles using the Data class, there is a flag to specify if that is
total reduced cross sectin or the charm contribution

Define Minuit parameters

Create DISFitter class, and pass Minuit parameters and Data classes to it

Ask MINUIT to minimize and do hard work.


Some remarks of the code

* All dimensionful quantities are expressed in units of GeV^n

* DISFitter (dis.cpp/hpp) is the actual class calculating photon-proton cross sections and
reduced cross sections, and eventually, chi^2. Integration over dipole size r is performed
there. Accuracy of this integral is controlled by constants at the beginning of file
dis.cpp

* IPsat dipole amplitude is described in ipsat.cpp
I have tried to design it such that one can later on add support for fluctuating proton
structure.

* DGLAP code is in the fortran files alphaS.f and LO_evolution_routine.f
There is also a C++ implementation of the same code which is used by default,
Note that in the Fortran version
LO_evol(X, Q2, gluon, coupling, Ag, lambdag)
actually now returns alphas(Q^2) * LO_evol
This way it remains easy to keep definition of alphas consistent in the Fortran and
C++ (IPsat class) parts

* If DGLAP solver is thread safe, chi^2 can be computed in parallel. Parallerization
is enabled by #define PARALLEL_CHISQR in dis.hpp

* Photon wave function is defined in virtual_photon.{cpp,hpp}.
In virtual_photon.hpp one can define #define USE_INTERPOLATOR,
which makes the code to calculate 2D grid of the wave function integrated
over z, and then to interpolate in that grid when calculating the cross sections.
This is not as accurate as numerically calulating the z integral at given r,Q^2
And may not be much faster, either...
Accuracy of the numerical z integration is controlled by constant defined at the beginning
of virtual_photon.cpp



