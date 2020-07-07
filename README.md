# Saras - Finite difference solver

Saras is an OpenMP-MPI hybrid parallelized Navier-Stokes equation solver written in C++.
It uses the finite-difference method for calculating spatial derivatives and parallelized geometric multigrid method for solving
the pressure Poisson equation.

All the source and library files for the Saras solver are contained in the following directories:

* ``./src/`` - contains all the different solvers available in Saras
* ``./lib/`` - contains all the libraries used by the solvers in ``./src/`` folder
* ``./compile/`` - contains the installation script to build and compile the solver
* ``./output/`` - the solution files written by the solver are contained in this folder, it also contains Python post-processing scripts

## Installing SARAS

To install ``SARAS``, you need to first clone the git repository into your local machine

`git clone https://github.com/roshansamuel/saras.git`

On LINUX systems which use the Bash shell, ``SARAS`` can be compiled by simply running the ``compileSaras.sh`` shell script in the `compile/` folder, as below

`bash compileSaras.sh`

The first few lines of the ``compileSaras.sh`` script can be used to set certain compilation parameters and flags:

* ``PROC`` - Number of processors to be used when running ``SARAS``. This parameter is used only if the ``EXECUTE_AFTER_COMPILE`` parameter is uncommented.
* ``REAL_TYPE`` - ``SARAS`` supports computations with both double and single precision floating point values. This parameter must be either ``SINGLE`` or ``DOUBLE``
* ``PLANAR`` - This parameter has to be enabled to use the ``SARAS`` executable for 2D simulations.
* ``TIME_RUN`` - Suppresses file-writing and I/O operations. This flag is enabled only when timing the solver for scaling runs.
* ``TEST_RUN`` - Runs the unit-tests module of the solver. This is distinct from the tests contained in the ``tests/`` folder of the solver.
* ``EXECUTE_AFTER_COMPILE`` - The script automatically runs the executable by issuing the ``mpirun`` command. This flag is enabled mainly during development for quickly compiling and running the solver.

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

Currently, ``SARAS`` is compatible only with yaml-cpp 0.3 and below, since it uses ``YAML::Parser::GetNextDocument``,
which was removed in subsequent versions.
On Debian based systems, the YAML library can be installed by running

`sudo apt-get install libyaml-cpp-dev`

However, this may install yaml-cpp 0.5 or newer, and the older version (which is still in available in many repositories)
has to be specifically installed for ``SARAS`` to run.
Otherwise, the compatible version of yaml-cpp library can be downloaded from [here](http://turbulencehub.org/wp-content/uploads/Download_Files/yaml-cpp-release-0.3.0.tar.gz).
Please extract the archive and follow the installation instructions.

### HDF5

The HDF5 library has to be installed with parallel writing enabled for the file writing functions of ``SARAS`` to work.
The library can either be downloaded and manually installed from [here](http://turbulencehub.org/wp-content/uploads/Download_Files/hdf5-1.8.20.tar.bz2),
or the native package manager of the OS can be used to locate and install the library.

More instructions on installing the libraries listed above can be found [here](http://turbulencehub.org/index.php/codes/tarang/installing-tarang/).
If any of the above libraries is being installed to the home directory,
please make sure to update the relevant paths in the shell configuration file (``~/.bashrc``)

## Running SARAS

``SARAS`` can be executed by issuing the ``mpirun`` command at the root folder of the solver (assuming that MPICH is installed as mentioned above).

``mpirun -np <number_of_processors> ./saras``

It is essential to set the parameters appropriately with the ``parameters.yaml`` file in the ``input/`` folder of the solver.
The number of processors specified to the ``mpirun`` command should be equal to the product of ``X Number of Procs`` and ``Y Number of Procs`` options
within the ``Parallel`` sub-section of ``parameters.yaml``.
Please check the ``parameters.yaml`` file for the full list of options specifiable to the solver, and their explanations (in comments).

For more information please refer to the ``SARAS`` [documentation](https://roshansamuel.github.io/saras/).

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

At the end of the test, a plot of the x and y velocity profiles is shown to the user and saved as ``ldc_validation.png`` in the folder ``tests/ldcTest/``.

## Setting up a new case in SARAS

2D and 3D cases require separate executables of ``SARAS``.
Once the dimensionality of the case is decided, ``SARAS`` has to compiled with appropriate compilation flags.
The executable file may then be shifted to a work folder where it will be run.

### Get the executable file

Users on UNIX systems are encouraged to use the BASH script, ``compileSaras.sh``, in the ``./compile/`` folder to compile ``SARAS``.
Immediately below the license header at the top of the script, there are 6 flags which the user can enable/disable.
A flag can be enabled or disabled by uncommenting or commenting it respectively.

It is best to disable the ``EXECUTE_AFTER_COMPILE`` flag so that script doesn't run the executable immediately after compilation.
If the ``EXECUTE_AFTER_COMPILE`` flag is disabled, the ``PROC`` variable, which sets the number of processes to run the executable with, can be ignored.
The user must enable or disable the ``PLANAR`` flag depending on whether the case is 2D or 3D respectively.

The ``REAL_TYPE`` flag has the default value of ``DOUBLE``.
This indicates that the executable will use double precision floating point numbers.
To use single precision floats, set the flag to ``SINGLE``.
The remaining two flags, ``TIME_RUN`` and ``TEST_RUN``, are used in special cases only and is best left to their default disabled states.

Once the above shell script variables have been set, the script can be executed at the command line to compile ``SARAS``.
If the compilation occurs without hiccups, an executable file named ``saras`` will appear in the root folder of the solver.

### Set the parameters file

When ``saras`` is executed, it will first read the case parameters from a YAML file named ``parameters.yaml``.
The user must set these parameters appropriately before executing ``saras``.

A sample ``parameters.yaml`` file is provided with the solver in the ``./input/`` folder.
The parameters are grouped under 5 sections, viz., ``Program``, ``Mesh``, ``Parallel``, ``Solver`` and ``Multigrid``.

* The main parameters that set the boundary conditions, initial conditions, forcing/source terms, etc. are found under the ``Program`` section.
* Grid parameters like number of points, stretching parameter for non-uniform grids, etc. are found under the ``Mesh`` section.
* ``Parallel`` section lets the user define how many MPI sub-domains to decompose the computational domain into, and the number of OpenMP threads to use.
* Non-dimensional time-step, file write intervals, final non-dimensional time and so on are set in the ``Solver`` section.
* Finally, ``Multigrid`` section lets the user tweak the parameters of the Geometric Multigrid solver used to solve the pressure Poisson equation.

Each parameter has documentation written into the ``parameters.yaml`` file itself.

> Additionally, the user can add custom boundary conditions to the ``boundary`` library.
> The source files of the ``boundary`` library can be found in ``./lib/boundary`` folder of the solver.
> Similarly, the user can add custom initial conditions to the ``initial`` library in the ``./lib/initial/`` folder,
> and custom forcing/source terms to the ``force`` library in ``./lib/force/``.
> All the source files in these libraries have extensive Doxygen documentation, and are written to be as self-explanatory as possible.

### Running and processing data

The folder where the executable ``saras`` will be run must contain two subfolders - ``./input/`` and ``./output/``.
The ``parameters.yaml`` must be saved in ``./input/``, and the solver will write data into ``./output/``.

Based on the values in ``parameters.yaml``, the solver will write solution data, time series, probe measurements, etc. in the ``./output/`` folder.
The solver will also periodically dump the entire field data into a file named ``restartFile.h5``, in the ``./output/`` folder.
This file will be read by the solver to resume computations, should it stop before completing the simulation.

The solution data is written in HDF5 format, while time-series and probe data are written in ASCII format.
Many open source visualization software are capable of reading HDF5 data format.
Moreover, Python can also read HDF5 files using the ``h5py`` module.

## License

``SARAS`` is an open-source package made available under the New BSD License.

## Contributions and bug reports

Contributions to this project are very welcome.
If you wish to contribute, please create a branch with a [pull request](https://github.com/roshansamuel/saras/pulls) and the proposed changes can be discussed there.

If you find a bug, please open a new [issue](https://github.com/roshansamuel/saras/issues/new) on the GitHub repository to report the bug.
Please provide sufficient information for the bug to be reproduced.

## References

Various articles and pages used to make programming decisions during development of the solver are listed here:

### General articles

1. https://stackoverflow.com/questions/4816698/avoiding-circular-dependencies-of-header-files
2. https://stackoverflow.com/questions/8111677/what-is-argument-dependent-lookup-aka-adl-or-koenig-lookup
3. https://www.codesynthesis.com/~boris/blog/2012/04/04/when-provide-empty-destructor/

### Articles on multi-grid methods

1. http://math.mit.edu/classes/18.086/2006/am63.pdf
2. http://www.mgnet.org/mgnet-tuts.html

### Journal references

1. Ghia, U., Ghia, K. N., & Shin, C. T. (1982). High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method. J. Comput. Phys., 48(3), 387-411. 
