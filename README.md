# Saras - Finite difference solver

Saras is an OpenMP-MPI hybrid parallelized Navier-Stokes equation solver written in C++.
It uses the finite-difference method for calculating spatial derivatives and parallelized geometric multi-grid method for solving
the pressure Poisson equation.

All the source and library files for the Saras solver are contained in the following directories:

* ``./src/`` - contains all the different solvers available in Saras
* ``./lib/`` - contains all the libraries used by the solvers in ``./src/`` folder
* ``./compile/`` - contains the installation script to build and compile the solver
* ``./output/`` - the solution files written by the solver are contained in this folder, it also contains Python post-processing scripts

## Installing SARAS

``SARAS`` relies on a few libraries for its calculations.
Therefore the first step towards building ``SARAS`` is to install the following dependencies:

* ``blitz`` - All array manipulations are performed using the Blitz++ library
* ``cmake`` - Necessary to build the ``saras`` executable
* ``mpich`` - For parallel computation using MPI
* ``yaml`` - The input parameters are stored in a YAML file which needs to be parsed using yaml-cpp library.
* ``hdf5`` - The output files are written in HDF5 format

Packages like ``cmake``, ``mpich`` and ``hdf5`` can be installed from the OS package manager,
while ``yaml`` and ``blitz`` can be downloaded and installed manually.

However, we provide below an alternative series of steps that will install all the libraries in the user's home directory.
This method has three advantages:

* It does not require administrator (``sudo``) privileges,
* It will not disturb pre-existing packages already installed on the system
* The method has been tested on both Linux (Ubuntu 14.04 and above), and MacOS (Mojave)

### Download all the dependencies

All the required packages can be downloaded from our [lab website](http://turbulencehub.org), and installed over a terminal.
It is advisable to create a temporary directory where the packages can be downloaded and extracted.
After navigating to the temporary directory, download the packages using ``wget`` on Linux, or ``curl`` on MacOS.
For instance, on MacOS, the following lines will download ``CMake``, ``Blitz++``, ``yaml-cpp``, ``MPICH`` and ``HDF5`` packages respectively:

```
curl -O https://turbulencehub.org/wp-content/uploads/Download_Files/cmake-2.8.12.tar.gz
curl -O https://turbulencehub.org/wp-content/uploads/Download_Files/blitz-1.0.1.tar.gz
curl -O https://turbulencehub.org/wp-content/uploads/Download_Files/yaml-cpp-release-0.3.0.tar.gz
curl -O https://turbulencehub.org/wp-content/uploads/Download_Files/mpich-3.1.3.tar.gz
curl -O https://turbulencehub.org/wp-content/uploads/Download_Files/hdf5-1.8.20.tar.bz2
```

On Linux, please replace ``curl -O`` with ``wget``.
If you are installing on a remote machine which doesn't have direct internet access,
the packages can be downloaded using the links listed above, and transferred over ssh or ftp.

### Extract all the packages

The ``tar`` command (available by default on both Linux and MacOS) can be used to extract the above archives.
For instance, ``CMake`` can be extracted by the command:

`tar -xf cmake-2.8.12.tar.gz`

On MacOS, if the above command fails with the error ``Failed to set default locale``, the $LANG environment variable has to be set as shown below:

`export LANG=en_US.UTF-8`

Once all the downloaded packages have been extracted, create an install location in the user's home folder if it doesn't exist already.
Usually, packages are installed in ``$HOME/local/`` directory.
The next steps will assume that the folder ``local/`` exists in the user's home directory.

### Install CMake

Navigate to the folder created by extracting the ``CMake`` package, configure the installation script, and install the package:

```
cd cmake-2.8.12/
./configure --prefix=$HOME/local
make -j4 install
```

Note that the option ``-j4`` given to ``make`` will use 4 cores of your system.
If more cores are available, the process can be speeded-up by specifying a higher number.

### Export path variables

Once ``CMake`` has been installed in ``$HOME/local`` it should be made available for the next steps of installation.
For this, the path variables have to be updated to let the system know that ``cmake`` is installed in ``$HOME/local``.
Accordingly, set the following environment variables:

```
export PATH=$HOME/local/bin:$PATH
export PKG_CONFIG_PATH=$HOME/local/lib/pkgconfig:$PKG_CONFIG_PATH
export HDF5_ROOT=$HOME/local
export CPATH=$HOME/local/include/:$CPATH
export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=$HOME/local/lib:$LIBRARY_PATH
export MANPATH=$HOME/local/share/man/:$MANPATH
```

Please note that the environment variables set here (and the $LANG variable set above if needed), exist only for the duration of the terminal session.
If the session is terminated by closing the terminal or logging out, the variables will be reset.
To permanently add the path variables, the above lines may be appended to the shell profile file (like ``bashrc`` or ``bash_profile``).

### Install Blitz++

Navigate to the folder created by extracting the ``Blitz++`` package, and install the package as done for ``CMake``:

```
cd ../blitz-1.0.1/
./configure --prefix=$HOME/local
make -j4 install
```

### Install yaml-cpp

Installing the ``yaml-cpp`` package will require ``cmake``, and hence uses a slightly different syntax.

```
cd ../yaml-cpp-release-0.3.0/
cmake -DCMAKE_INSTALL_PREFIX=$HOME/local
make -j4 install
```

With ``cmake``, the ``-DCMAKE_INSTALL_PREFIX`` argument performs the same function as ``--prefix`` for ``make``, namely specifying the install directory.

Currently, ``SARAS`` is compatible only with yaml-cpp 0.3 and below, since it uses ``YAML::Parser::GetNextDocument``, which was removed in subsequent versions.
If ``yaml-cpp`` is installed using a package manager, the OS may install yaml-cpp 0.5 or newer.
The older version (which is still in available in many repositories) has to be specifically installed for ``SARAS`` to run.

### Install MPICH

Repeat the same steps used to install ``cmake`` and ``blitz`` from the folder into which the ``MPICH`` package was extracted:

```
cd ../mpich-3.1.3/
./configure --prefix=$HOME/local
make -j4 install
```

### Install HDF5 library

The ``HDF5`` library requires ``MPICH`` so that it can perform parallel file I/O operations.
Hence it must be installed only after installing ``MPICH`` as done above.
A few additional build flags are also provided when building the ``HDF5`` library:

```
cd ../hdf5-1.8.20/
CC=mpicc CXX=mpicxx ./configure --prefix=$HOME/local --enable-parallel --without-zlib
make -j4 install
```

Note that while building the ``HDF5`` library, it is being explicitly specified that the MPI compiler must be used.

### Clone and install SARAS

With luck, the above steps will have installed all the dependencies required by ``SARAS``.
They have been tested with ``gcc`` on Linux, ``homebrew-gcc`` as well as ``clang`` on MacOS.
Now the ``saras`` repository can be cloned into your machine, and compiled:

```
git clone https://github.com/roshansamuel/saras.git
cd saras/compile/
bash compileSaras.sh
```

The build script, ``compileSaras.sh`` automatically builds ``SARAS`` with a few default parameters.
It has been tested to work on both Linux and MacOS.
If the compilation is successful, the ``saras`` executable will be built and linked in the root ``saras/`` folder.

Note that the first few lines of the ``compileSaras.sh`` script can be used to set certain compilation parameters and flags:

* ``PROC`` - Number of processors to be used when running ``SARAS``. This parameter is used only if the ``EXECUTE_AFTER_COMPILE`` parameter is uncommented.
* ``REAL_TYPE`` - ``SARAS`` supports computations with both double and single precision floating point values. This parameter must be either ``SINGLE`` or ``DOUBLE``
* ``PLANAR`` - This parameter has to be enabled to use the ``SARAS`` executable for 2D simulations.
* ``TIME_RUN`` - Suppresses file-writing and I/O operations. This flag is enabled only when timing the solver for scaling runs.
* ``EXECUTE_AFTER_COMPILE`` - The script automatically runs the executable by issuing the ``mpirun`` command. This flag is enabled mainly during development for quickly compiling and running the solver.


## Running SARAS

``SARAS`` can be executed by issuing the ``mpirun`` command at the root folder of the solver.

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
Additionally, the convergence of the multi-grid Poisson solver of ``SARAS`` can also be tested.
This test is also available in the ``tests/`` folder.

## Setting up a new case in SARAS

2D and 3D cases require separate executables of ``SARAS``.
The executable has to be placed in a folder with two sub-folders: ``input/`` and ``output/``.
``SARAS`` will read parameters from the ``input/`` folder, and write solution data into the ``output/`` folder.
Following are the steps to get a case up and running:

### Get the executable file

Based on the dimensionality of the problem being solved, and the precision of floating point numbers to be used,
set the ``PLANAR`` and ``REAL_TYPE`` variables in the build script - ``compileSaras.sh``.
Executing this shell script will produce the ``saras`` executable file as mentioned above.
It is best to disable the ``EXECUTE_AFTER_COMPILE`` flag so that the script doesn't run the executable immediately after compilation.

### Set the parameters file

When ``saras`` is executed, it will first read the case parameters from a YAML file named ``parameters.yaml``.
The user must set these parameters appropriately before executing ``saras``.

A sample ``parameters.yaml`` file is provided with the solver in the ``./input/`` folder.
The parameters are grouped under 5 sections, viz., ``Program``, ``Mesh``, ``Parallel``, ``Solver`` and ``Multigrid``.

* The main parameters that set the boundary conditions, initial conditions, forcing/source terms, etc. are found under the ``Program`` section.
* Grid parameters like number of points, stretching parameter for non-uniform grids, etc. are found under the ``Mesh`` section.
* ``Parallel`` section lets the user define how many MPI sub-domains to decompose the computational domain into, and the number of OpenMP threads to use.
* Non-dimensional time-step, file write intervals, final non-dimensional time and so on are set in the ``Solver`` section.
* Finally, ``Multigrid`` section lets the user tweak the parameters of the Geometric Multi-grid solver used to solve the pressure Poisson equation.

Each parameter has documentation written into the ``parameters.yaml`` file itself.

> Additionally, the user can add custom boundary conditions to the ``boundary`` library.
> The source files of the ``boundary`` library can be found in ``./lib/boundary`` folder of the solver.
> Similarly, the user can add custom initial conditions to the ``initial`` library in the ``./lib/initial/`` folder,
> and custom forcing/source terms to the ``force`` library in ``./lib/force/``.
> All the source files in these libraries have extensive Doxygen documentation, and are written to be as self-explanatory as possible.

### Running and processing data

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
