# Saras - Finite difference solver

Saras is an OpenMP-MPI hybrid parallelized Navier-Stokes equation solver written in C++.
It uses the finite-difference method for calculating spatial derivatives and parallelized algebraic multi-grid method for solving
the pressure Poisson equation.

All the source and library files for the Saras solver are contained in the following directories:

* ``./src/`` - contains all the different solvers available in Saras
* ``./lib/`` - contains all the libraries used by the solvers in ``./src/`` folder
* ``./compile/`` - contains the installation script to build and compile the solver
* ``./output/`` - the solution files written by the solver are contained in this folder, it also contains Python post-processing scripts

## Installation

To install ``SARAS``, you need to first clone the git repository into your local machine

`git clone https://github.com/roshansamuel/saras.git`

On LINUX systems which use the Bash shell, ``SARAS`` can be compiled by simply running the ``compileSaras.sh`` shell script in the `compile/` folder, as below

`bash compileSaras.sh`

Before compilation, a few dependencies have to installed first.

* ``blitz`` - All array manipulations are performed using the Blitz++ library
* ``cmake`` - Necessary to make the executable from the source files
* ``mpich`` - The compiler used is mpic++
* ``yaml`` - The input parameters are stored in the parameters.yaml file which needs the yaml-cpp library to parse.
* ``hdf5`` - The output files are written in HDF5 format

### Blitz++

To install the Blitz++ library, please download the library from [here](http://turbulencehub.org/wp-content/uploads/Download_Files/blitz-1.0.1.tar.gz), and extract the archive.
Follow the installation instructions in the archive.

### CMake

On LINUX systems which use Debian package manager, please use the package manager itself to install CMake.
For example, on Ubuntu systems, CMake can be installed by running

`sudo apt-get install cmake`

### MPICH

Similar to CMake installation above, it is best to install MPICH using the native package manager.

### YAML

On Debian based systems, the YAML library can be installed by running

`sudo apt-get install libyaml-cpp-dev`

Otherwise, the yaml-cpp library can be downloaded from [here](http://turbulencehub.org/wp-content/uploads/Download_Files/yaml-cpp-release-0.3.0.tar.gz).
Please extract the archive and follow the installation instructions.

### HDF5

The HDF5 library has to be installed with parallel writing enabled for the file writing functions of ``SARAS`` to work.
The library can either be downloaded and manually installed from [here](http://turbulencehub.org/wp-content/uploads/Download_Files/hdf5-1.8.20.tar.bz2),
or the native package manager of the OS can be used to locate and install the library.


More instructions on installing the libraries listed above can be found [here](http://turbulencehub.org/index.php/codes/tarang/installing-tarang/).
If any of the above libraries is being installed to the home directory,
please make sure to update the relevant paths in the shell configuration file (``~/.bashrc``)

## Testing SARAS

``SARAS`` offers an automated testing process to validate the solver after installation.
The relevant test scripts can be found in the ``tests/`` folder of the solver.
Executing the Bash shell script ``testSaras.sh``, will compile ``SARAS``, and run it with a pre-defined set of parameters.
We use the benchmark results on 2D lid-driven cavity (LDC) performed by Ghia et al (1982) to validate ``SARAS``.
The test can be executed by running the following command within the ``tests/`` folder.

``bash testSaras.sh``

The test uses 4 cores and takes about 12 minutes to complete on an Intel workstation.
At the end of the test, the Python script ``validate_ldc.py``, found in ``tests/ldcTest/`` reads the output from ``SARAS``,
and plots the velocity profiles along with the data from Ghia et al's result.

The following Python modules are necessary for the Python test script to execute successfully

* numpy
* matplotlib
* h5py
* yaml

At the end of the test, a plot of the x and y velocity profiles is drawn into a ``test.png`` file for the user to verify.
A measure of the deviation of the results from the expected values is also calculated and the test is deemed passed if this
deviation is within the set tolerance.

## License

``SARAS`` is an open-source package made available under the New BSD License.

## References

Various articles and pages used to make programming decisions during development of the solver are listed here:

### General stuff

1. http://coding.derkeiler.com/Archive/C_CPP/comp.lang.c/2004-02/1382.html
2. https://stackoverflow.com/questions/4816698/avoiding-circular-dependencies-of-header-files
3. http://blitz.sourceforge.net/resources/blitz-0.9.pdf
4. https://stackoverflow.com/questions/8111677/what-is-argument-dependent-lookup-aka-adl-or-koenig-lookup
5. https://www.codesynthesis.com/~boris/blog/2012/04/04/when-provide-empty-destructor/

### Articles on multi-grid methods

1. http://math.mit.edu/classes/18.086/2006/am63.pdf
2. http://www.mgnet.org/mgnet-tuts.html

### Journal references

1. Ghia, U., Ghia, K. N., & Shin, C. T. (1982). High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method. J. Comput. Phys., 48(3), 387-411. 
