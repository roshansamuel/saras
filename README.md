# Saras - Finite difference solver

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02095/status.svg)](https://doi.org/10.21105/joss.02095)

SARAS is an MPI parallelized Navier-Stokes equation solver written in C++.
It uses the finite-difference method for calculating spatial derivatives and parallelized geometric multi-grid method for solving
the pressure Poisson equation.

All the relevant files of the solver are contained in the following directories:

* ``./src/`` - contains the different solvers available in SARAS
* ``./lib/`` - contains all the libraries used by the solvers in ``./src/`` folder
* ``./input/`` - contains the parameters to be read by the solver
* ``./output/`` - the solution files are written into this folder, it also contains Python post-processing scripts

## Installing SARAS

``SARAS`` relies on a few libraries for its calculations.
Therefore the first step towards building ``SARAS`` is to install the following dependencies:

* ``cmake`` - Necessary to build the ``saras`` executable
* ``mpich`` - For parallel computation using MPI
* ``hdf5`` - The output files are written in HDF5 format
* ``blitz`` - All array manipulations are performed using the Blitz++ library
* ``yaml-cpp`` - The input parameters are stored in a YAML file which needs to be parsed using yaml-cpp library.

Packages like ``cmake``, ``mpich``, and ``hdf5`` can be installed from the OS package manager.
However, the latest versions of ``blitz`` and ``yaml-cpp`` packages can be downloaded from GitHub as described below.

If you do not have sudo privileges to install packages using the OS package manager,
you can install them in your home folder instead.
This also offers the potential advantage of not disturbing pre-existing packages already installed on the system.
The steps listed below explain this method of installation.

That said, libraries like ``cmake``, ``MPICH`` and ``HDF5`` are normally available on most computing systems.
Please check if these packages are already installed, and if they are, you can skip their installation steps.

### Download all the dependencies

All the required packages can be downloaded and installed over a terminal.
It is advisable to create a temporary directory where the packages can be downloaded and extracted.
After navigating to the temporary directory, download the packages using ``wget`` command on Linux, or ``curl`` command on MacOS.
[CMake](https://cmake.org/download/),
[MPICH](https://www.mpich.org/downloads/) and
[HDF5](https://portal.hdfgroup.org/display/support/HDF5+1.8.21)
packages can be downloaded from their respective sites, and extracted by the ``tar`` command.

```
wget https://cmake.org/files/v3.12/cmake-3.12.4.tar.gz
wget http://www.mpich.org/static/downloads/3.4.2/mpich-3.4.2.tar.gz
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.21/src/hdf5-1.8.21.tar.gz

tar -xf cmake-3.12.4.tar.gz
tar -xf mpich-3.4.2.tar.gz
tar -xf hdf5-1.8.21.tar.gz
```

On Mac OS, please replace ``wget`` with ``curl -O``.
It is best to download the latest versions of ``Blitz++`` and ``yaml-cpp`` from their respective Git repositories.

```
git clone https://github.com/jbeder/yaml-cpp.git
git clone https://github.com/blitzpp/blitz.git
```

If you are installing on a remote machine which doesn't have direct internet access,
the packages can be downloaded using the links listed above, and transferred over ssh or ftp.

> On MacOS, if the ``tar`` command fails with the error ``Failed to set default locale``,
> the $LANG environment variable has to be set.
> For this you need to execute the command: `export LANG=en_US.UTF-8`

Once all the downloaded packages have been extracted, create an install location in your home folder if it doesn't exist already.
Usually, packages are installed in ``$HOME/local/`` directory.
The next steps will assume that the folder ``local/`` exists in the user's home directory.

### Install CMake

After navigating to the folder created by extracting the ``CMake`` package,
execute the following commands to configure the installation script and install the package:

```
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
export CPATH=$HOME/local/include/:$CPATH
export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=$HOME/local/lib:$LIBRARY_PATH
export MANPATH=$HOME/local/share/man/:$MANPATH
```

Please note that the environment variables set here exist only for the duration of the terminal session.
If the session is terminated by closing the terminal or logging out, the variables will be reset.
To permanently add the path variables, the above lines may be appended to the shell profile file (like ``bashrc`` or ``bash_profile``).

### Install yaml-cpp

Installing the ``yaml-cpp`` package will require ``cmake``.
It is best to build the package within a temporary ``build/`` directory to avoid disturbing the source files.
After navigating to the folder created by cloning the ``yaml-cpp`` repository, execute the following steps:

```
mkdir build
cd build/
cmake -DCMAKE_INSTALL_PREFIX=$HOME/local ../
make -j4 install
```

With ``cmake``, the ``-DCMAKE_INSTALL_PREFIX`` argument performs the same function as ``--prefix`` for ``make``, namely specifying the install directory.

> ``SARAS`` is now compatible with the latest versions of ``yaml-cpp``.
> The older ``yaml-cpp 0.3`` uses ``YAML::Parser::GetNextDocument`` to parse the YAML file.
> However, this function is not secure and hence deprecated in later versions.
> If you have an older version of ``yaml-cpp`` installed on your system, you can still run ``SARAS`` with the older package.
> To enable this, the user has to pass the ``-DYAML_LEGACY`` flag to ``CMake`` when building ``SARAS`` later.

### Install MPICH

Similar to the installation of ``CMake``, ``MPICH`` also uses a configure script to install the package.
After navigating to the folder created by extracting the ``mpich`` tarball, configure and install the package as done before:

```
./configure --prefix=$HOME/local
make -j4 install
```

> If the ``configure`` step throws the error ``The Fortran compiler gfortran will not compile files that call the same routine with arguments of different types.``,
> please set the following flag, and rerun the configure script.
>
> `export FFLAGS="-w -fallow-argument-mismatch -O2"`

### Install Blitz++

Although older versions of ``Blitz`` used a configure script to build the package, the latest versions use ``CMake`` instead.
Similar to the installation of ``yaml-cpp``, you will have to build the package within a temporary ``build/`` directory.
After navigating to the folder created by cloning the ``blitz`` repository, execute the following commands:

```
mkdir build
cd build/
cmake -DCMAKE_INSTALL_PREFIX=$HOME/local ../
make -j4 install
```

The latest version of ``CMake`` is required to build the ``blitz`` package.
If the above steps fail due to ``CMake`` version mismatch, you will have to upgrade your ``CMake``.
If you installed ``CMake`` as done in the steps above, you will not face this issue,
since the latest version of ``CMake`` was downloaded for the installation process described above.

> At the final step in ``make install``, the Blitz++ installer may require ``python2``.
> If you encounter an error due to Python version mismatch,
> please switch the environment to use ``python2`` temporarily for this step of the installation.

### Install HDF5 library

The ``HDF5`` library requires ``MPICH`` so that it can perform parallel file I/O operations.
Hence it must be installed only after installing ``MPICH`` (if ``MPICH`` was not already available on your system).
A few additional build flags are also provided when building the ``HDF5`` library:

```
CC=mpicc CXX=mpicxx ./configure --prefix=$HOME/local --enable-parallel --without-zlib
make -j4 install
```

> Note that while building the ``HDF5`` library, it is being explicitly specified that the MPI compiler must be used.
> On some systems, it might be necessary to add an extra compiler flag before running ``make`` in order to install ``HDF5`` properly:
> 
> `export CFLAGS=-Wno-error=implicit-function-declaration`

### Clone and install SARAS

With luck, the above steps will have installed all the dependencies required by ``SARAS``.
Now the ``saras`` repository can be cloned into your machine:

```
git clone https://github.com/roshansamuel/saras.git
```

After navigating to the folder created by cloning the repository,
you can follow the same steps for building the package as done when installing ``Blitz``, ``yaml-cpp``, etc.

```
mkdir build
cd build/
CC=mpicc CXX=mpicxx cmake ../
make -j4
```

``CMake`` will build the ``saras`` executable in the root folder of the solver.
The following flags can be passed to ``CMake`` to enable/disable different features of ``SARAS``:

* ``-DPLANAR=ON`` - This compiles ``SARAS`` for 2D simulations. By default, ``SARAS`` is compiled for 3D runs
* ``-DREAL_SINGLE=ON`` - ``SARAS`` can compute in single-precision. By default, ``SARAS`` uses double-precision
* ``-DTIME_RUN=ON`` - This flag suppresses file-writing and I/O when compiling the solver for scaling studies
* ``-DYAML_LEGACY=ON`` - As described previously, this flag allows ``SARAS`` to use older versions of ``yaml-cpp``

For example, to build ``SARAS`` for a 2D simulation using single-precision calculations, the solver will be configured as:

```
CC=mpicc CXX=mpicxx cmake -DPLANAR=ON -DREAL_SINGLE=ON ../
```

Note that the MPI compilers have to be specified to ``CMake`` for the MPI headers to be found when building the executable.

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
Executing the Bash shell script ``testLDC.sh``, will compile ``SARAS``, and run it with a pre-defined set of parameters.
We use the benchmark results on 2D lid-driven cavity (LDC) performed by Ghia et al (1982) to validate ``SARAS``.
The test can be executed by running the following command within the ``tests/`` folder.

``bash testLDC.sh``

The test uses 4 cores and takes about 12 minutes to complete on an Intel workstation.
At the end of the test, the Python script ``checkLDC.py``, found in ``tests/ldcTest/`` reads the output,
and plots the velocity profiles along with the data from Ghia et al's result.

The following Python modules are necessary for the Python test script to execute successfully

* numpy
* scipy
* matplotlib
* h5py
* yaml

At the end of the test, a plot of the x and y velocity profiles is shown to the user and saved as ``ldc_validation.png`` in the folder ``tests/ldcTest/``.
Additionally, the convergence of the multi-grid Poisson solver of ``SARAS`` can also be tested.
This test is also available in the ``tests/`` folder.

## Setting up a new case in SARAS

After generating the ``saras`` executable file as described in the section on [Installing SARAS](#clone-and-install-saras),
The executable has to be placed in a folder with two sub-folders: ``input/`` and ``output/``.
``SARAS`` will read parameters from the ``input/`` folder, and write solution data into the ``output/`` folder.

### Set the parameters file

The parameters of a case to be simulated with ``SARAS`` are specified in a YAML file named ``parameters.yaml``.
The user must set these parameters appropriately before executing ``saras``.

A sample ``parameters.yaml`` file is provided with the solver in the ``./input/`` folder.
The parameters are grouped under 5 sections, viz., ``Program``, ``Mesh``, ``Parallel``, ``Solver`` and ``Multigrid``.

* The main parameters that set the boundary conditions, initial conditions, forcing/source terms, etc. are found under the ``Program`` section.
* Grid parameters like number of points, stretching parameter for non-uniform grids, etc. are found under the ``Mesh`` section.
* ``Parallel`` section lets the user define how many MPI sub-domains to decompose the computational domain.
* Non-dimensional time-step, file write intervals, final non-dimensional time and so on are set in the ``Solver`` section.
* Finally, ``Multigrid`` section lets the user tweak the parameters of the Geometric Multi-grid solver used to solve the pressure Poisson equation.

Each parameter has documentation written into the ``parameters.yaml`` file itself.

> Additionally, the user can add custom boundary conditions to the ``boundary`` library.
> The source files of the ``boundary`` library can be found in ``./lib/boundary`` folder of the solver.
> Similarly, the user can add custom initial conditions to the ``initial`` library in the ``./lib/initial/`` folder,
> and custom forcing/source terms to the ``force`` library in ``./lib/force/``.
> All the source files in these libraries have extensive Doxygen documentation, and are written to be as self-explanatory as possible.

### Running and processing data

The solver will write the solution data files into ``./output/`` folder.
Based on the values in ``parameters.yaml``, the solver may write solution data, time series, probe measurements, etc. in this folder.
The solver will also periodically dump the entire field data into a file named ``restartFile.h5``.
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
