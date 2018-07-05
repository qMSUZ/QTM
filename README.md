# The programming library implementing Quantum Trajectories Method => QTM package

The QTM package is a programming library implementing Quantum Trajectories Method with use of MPI protocol. The MPI standard allows to employ parallel computing what increases performance of calculation -- this is especially important for simulation of quantum systems, which complexity is exponential and many trajectories have to be calculated. 

The QTM package was prepared in the C++ language but also utilizes numerical methods, to solve Ordinary Differential Equations (ODEs), from the ZVODE package.

# Installation

The installation of QTM package using source code requires previously installed C++ and Fortran compilers. Also libraries Blas and Lapack must be installed. The last library is a package implementing MPI protocol, e.g. MPICH or OpenMPI.

The QTM package was developed within Ubuntu operating system (http://www.ubuntu.com), therefore below we present how to compile the package on Ubuntu.

The first step is an installation of basic tool packages for code developers:

sudo apt-get install build-essential

The next step is installing BLAS and LAPACK libraries:

sudo apt-get install libblas-dev liblapack-dev

Now a package for MPI, here OpenMPI:

sudo apt-get install libopenmpi-dev

The current version of QTM package may be downloaded from the repository:

git clone https://github.com/qMSUZ/QTM

Let us enter the main catalog:

cd QTM

Before the compilation of QTM package we have to install the ZVODE package (it is not attached to the QTM repository):

make zvode_download

The command to compile the example for trilinear Hamiltonian is:

make triham-ex

Running the package does not require the cluster of workstations based on MPI standard. The QTM package may be also used on a multi-core personal computer. So irrespectively the hardware, we call:

mpirun -n 5 triham-ex

what results with running five MPI processes -- one of them is a master node and other work as computational nodes.

# Basic example

The QTM package does not  work in an interactive mode. At the beginning of calculation it needs: the Hamiltonian, collapse operators, the expectation value and pointing the method (from the ZVODE package) to solve ODEs. However, these are only parameters that have to be supplied by the user who wants to solve a problem with use of QTM.

The package contains some examples and here we present one of them which refers to a unitary Hamiltonian. The definitions of basic data structures are:

#include "complexnum.h"
#include "qtm.h" 

const size_t COLLAPSE_OPERATORS = 1;
const size_t Ntrj = 1;
const size_t N = 100;
const size_t WAVEVECTOR_LEAD_DIM = 2;
const size_t WAVEVECTOR_LEAD_DIM_SQR = 4;

uMatrix< simpleComplex<double>, WAVEVECTOR_LEAD_DIM > c_ops[ COLLAPSE_OPERATORS ];
uVector< simpleComplex<double>, WAVEVECTOR_LEAD_DIM_SQR > collapse_operator;
uVector< simpleComplex<double>, WAVEVECTOR_LEAD_DIM_SQR > expect_operator;
uVector< simpleComplex<double>, WAVEVECTOR_LEAD_DIM_SQR > H;
simpleComplex<double> alpha[WAVEVECTOR_LEAD_DIM];

#include "qtm.cc"

extra_options opt;

The QTM is a template, therefore at the beginning of the source file we have to add two header files which implement the basic data types and QTM itself. We also have to add mpi_main function which is placed in qtm.cc file.

Within the main function we may run the following code (msc is an abbreviation for function make_simpleComplex creating complex numbers):

int r = 0;
co[0] = msc( 0.0, 0.0 );  co[1] = msc( 0.05, 0.0 );
co[2] = msc( 0.05, 0.0 ); co[3] = msc( 0.0, 0.0 );

eo[0] = msc( 1.0, 0.0); eo[1] = msc( 0.0, 0.0);
eo[2] = msc( 0.0, 0.0); eo[3] = msc(-1.0, 0.0);

alpha[0] = msc( 1.0, 0.0);
alpha[1] = msc( 0.0, 0.0);
	
H[0] = msc( -0.00125, 0.0); 
H[1] = msc( 0.0, -0.62831853);
H[2] = msc( 0.0, -0.62831853); 
H[3] = msc( -0.00125, 0.0);
	
c_ops[0].rows=2; c_ops[0].cols=2;
c_ops[0].m = co; 

opt.type_output = OUTPUT_FILE;
opt.only_final_trj = 1;
opt.ode_method = METADAMS;
opt.tolerance = 1e-7;
opt.file_name = strdup("output-data.txt");
opt.fnc = &myfex_fnc_f1;
	
r = mpi_main<N, Ntrj, WV_LD, WV_LD_SQR, 1>(argc, argv, 1, 0, 10, 1, 1, opt);

The source code may be compiled directly:

g++ unitary-mc-ex.cc rgen_lfsr113.cc zvode.o zgesl.o zgefa.o zgbsl.o zgbfa.o -o unitary-ex -lblas -lgfortran -lopenmpi

Now, we can run the example (on a single PC or within a cluster of workstations):

mpirun -n 9 unitary-ex

This command runs 9 MPI processes (one master node and 8 computational nodes).



