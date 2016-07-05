# HipGISAXS: QUICK STARTER GUIDE #
Bug (what bug...?) reporting: Email: `asarje@lbl.gov`

## TABLE OF CONTENTS ##
  1. Licensing
  2. Supported Platforms
    1. Operating System Platforms
    2. System Hardware (Compute Environment)
    3. System Hardware (Processor Architectures)
  3. Software Pre-requisites
  4. HipGISAXS Directory Layout
  5. Building the Binary and Library
  6. Running the Binary
  7. Inputs
  8. Using the Library in an Application
  9. Appendix
    1. Special System Specific Build Instructions
      1. On a generic Linux (x86-64) system
      2. [NERSC](http://www.nersc.gov) Systems
        1. Edison (Cray XC30, Intel Ivy Bridge)
        2. Cori Phase 1 (Cray XC40, Intel Haswell)
      3. [OLCF](https://www.olcf.ornl.gov) Systems
        1. Titan (Cray XK7)
      4. For users on Andromeda/Bragg servers
    2. Installing Software Dependencies
      1. TIFF Library


## Licensing ##
The HipGISAXS software is only available to be downloaded and used by employees
of academic research institutions, not-for-profit research laboratories, or
governmental research facilities. Please read the accompanying LICENSE file
before downloading the software. By downloading the software, you are agreeing
to be bound by the terms of this Non-Commercial End User License Agreement.
This licensing information is also provided as a PDF file along with the documentation.


## Supported Platforms
HipGISAXS has been successfully tested on the following platforms.

### A. Operating System Platforms (Software Environment)
1. GNU/Linux x86\_64: Ubuntu, Red Hat Linux, SUSE Linux, Cray Linux Environment (XK7, XE6, XC30, XC40).
2. Darwin x86\_64: Mac OS X (El Capitan, Yosemite).
3. You could try HipGISAXS on any UNIX based OS: generally it should work.
4. Windows: Will probably be supported in future.

### B. System Hardware (Compute Environment)
1. Generic x86 laptop/desktop/server.
2. Clusters/Supercomputers based on x86 processors.
3. Generic x86 laptop/desktop/server equipped with Nvidia GPUs.
4. Clusters/Supercomputers equipped with Nvidia GPUs as accelerators on each node.

### C. System Hardware (Processor Architectures)
1. Intel/AMD processors (64-bit).
2. Nvidia GPUs with compute capability >= 2.0.
3. Intel MIC architecture: support under development.


## Software Pre-requisites
This software uses several third-party libraries, and they need to be installed and available in order to compile and run HipGISAXS.
For ease, if possible use the installations already available on your system, if any. Alternatively, download and install them yourself.
The following are the dependencies of this software:

### Required Software
1. **GNU C/C++ compilers**, version >= 4.7, OR  
   **Intel C/C++ compilers**, version >= 15 (with GNU compatibility >= 4.7).  
   To compile with GPU support enabled, it is recommended to use GNU compilers. 
2. **Scons**, version >= 2.0.  
   Scons can be obtained from: [http://scons.org](http://scons.org)  
   *NOTE: Scons is based on Python, so you need to have Python installed too.*
3. **Boost C++ Libraries**.  
   Boost can be obtained from: [http://www.boost.org](http://www.boost.org)  
4. **Tiff Library** (libtiff).  
   Tiff library can be obtained from: [http://remotesensing.org/libtiff](http://remotesensing.org/libtiff)  
   *NOTE: `libtiff` should be built with C++ support enabled. Also, disable `jpeg` support unless you have JPEG libraries installed. Brief instructions on building `libtiff` are given in the appendix.*

### Optional Software
1. **Nvidia CUDA toolkit** version >= 6.0.  
   CUDA can be obtained from [Nvidia CUDA downloads](http://developer.nvidia.com/cuda/cuda-downloads)  
   *NOTE: CUDA required to enable GPU support.*  
   *NOTE: CUDA is NOT required if compiling for CPU-only version.*
2. **OpenMPI** compiled with GNU or Intel compilers which ever you are using to compile HipGISAXS. Version > 1.4.4.  
   OpenMPI can be obtained from its [website](http://www.open-mpi.org/software)  
   Alternative implementations of MPI may also be used, such as MPICH and MVAPICH.  
   *NOTE: MPI library is required to enable support for multi-node systems.*  
   *NOTE: MPI library is NOT required if compiling for single node/server/desktop/laptop.*
3. **Parallel HDF5 library**.  
   HDF5 can be obtained from the [HDF5 group website](http://www.hdfgroup.org/downloads)  
   *NOTE: HDF5 library is NOT required if you do not plan to use input files in this format.*  
   *NOTE: HDF5 depends on the zlib and szip libraries.*  
   **zlib** can be obtained from: [http://www.zlib.net](http://www.zlib.net)  
   **szip** can be obtained from: [http://www.hdfgroup.org/doc\_resource/SZIP](http://www.hdfgroup.org/doc\_resource/SZIP)  


## HipGISAXS Directory Layout

    HipGISAXS
    |- README.md  : Duh!
    |- LICENSE    : This contains all the licensing information for HipGISAXS.
    |- SConstruct : Scons file for installation.
    |- SConscript : Scons file for installation.
    |- bin        : This contains HipGISAXS binaries generated by compilation.
    |- build      : A few Makefiles for various systems are provided in this.
    |- data       : This provides some sample input shape definition files.
    |- doc        : This will contain detailed documentation of HipGISAXS.
    |- examples   : This contains some example input files in HiG format.
    |- extras     : This contains some extra stuff such as syntax highlighting config for Vim.
    |- include    : All the source headers are within this directory. It contains subdirectories.
    |- lib        : The HipGISAXS library, libhipgisaxs.a, is generated here.
    |- man        : This has the man pages.
    |- obj        : All the object files generated during build are stored here.
    |- samples    : This contains some samples on how to use the HipGISAXS library.
    |- src        : This is the main source code directory. It contains many subdirectories.


## Building the Binary and Library
To build the HipGISAXS application binary and static library, use `scons`. Make sure you are passing paths to all the dependencies through `--extrapath=` option of scons:

    $ scons --extrapath=<path1>,<path2>,<etc>

To enable GPU support, use `--with-cuda` option:

    $ scons --extrapath=<path1>,<path2>,<etc> --with-gpu

To enable MPI support, use `--with-mpi` option:

    $ scons --extrapath=<path1>,<path2>,<etc> --with-gpu --with-mpi

or,

    $ scons --extrapath=<path1>,<path2>,<etc> --with-mpi

The generated binary, `hipgisaxs`, will be located in the `bin` directory.
The generated static library, `libhipgisaxs.a`, will be located in the `lib` directory.

... and you are done. Go enjoy HipGISAXS!

NOTE: See Appendix at the end of this file for more detailed and customized building information.

## Running the Binary
1. Make sure all the paths (data, output) and other variables are correctly set in the input HiG file.
2. Execute the binary with the input file as an argument. For example,  

   ```
     $ ./bin/hipgisaxs examples/01-pyramid.hig
   ```

## Inputs
The HipGISAXS binary takes as input a file in HiG format.
Please refer to detailed HipGISAXS documentation for details of the HiG format. Detailed information on the inputs is available at https://webhipgisaxs.lbl.gov  
A few sample input files are located in the directory `examples`, with extensions `.hig`.
Update the input files as needed.   
The main components in the input file to update are the following:

1. Shape name defines the input filename containing triangulated shape data.
   It should point to the correct location of the file (relative path):  
   ```
       ...
       shape = {
    	   ...
    	   name = "data/flexrod.obj"
       } ...
   ```
   This file needs to be either in HDF5 format or the OBJ format.  
   Some sample shape files are provided in the directory "data", with extensions ".hd5".

2. Output locatation needs to be defined as pathprefix and runname in the computation object:
   ```
      ...
      computation = {
    	  ...
    	  pathprefix = ".",
    	  runname = "myflexrod",
    	  ...
      } ...
   ```
   The `pathprefix` is a relative path to a directory.  
   The `runname` is appended with a timestamp, and a directory by this resulting name is created in the directory specified by `pathprefix`.  
   All output files are stored in this generated directory.  

3. Resolution alters the final image resolution, and also affects the run time:
   ```
      ...
      computation = {
    	  ...
    	  resolution = [ 1000 500 ],
    	  ...
      } ...
   ```
   Obviously, simulations with lower resolutions finish faster.  
   Modify the provided input file templates as needed.


## Using the Library in an Application
   Please refer to the examples provided in the `samples` directory.
   It contains simple code which uses the HipGISAXS library.
   The corresponsing Makefile is also provided as a reference.


## Appendix

### Special System Specific Build Instructions

#### A. On a generic Linux (x86-64) system
1. Ensure all above prerequisites are available on the system.
2. Ensure all system environment variables are set accordingly to include the prerequisites.
3. Once all the required software are available, use the `scons` command to build the binary. Example:  
   ```
      $ scons --extrapath=$PATHS_TO_SOFTWARE --with-mpi --with-cuda
   ```

4. On successful build, the binary will be generated in the `bin` directory, called `hipgisaxs`.

#### B. [NERSC](http://www.nersc.gov) Systems

##### Edison (Cray XC30, Intel Ivy Bridge)
1. All required software, except `libtiff` are available as modules on Edison. An example set of modules you can load are given in the `build/modules.edison` file. You could just `source` this file:  
   ```
      $ source build/modules.edison
   ```

2. You will need to install `libtiff`. Please refer to the required software section above.
3. An example build command is given in the `build/build-edison.sh`. If needed, make sure the paths are correctly set, including your installation of `libtiff`. Since Edison requires cross compilation for its compute nodes, make sure the `CC` and `CXX` environment variables are also set. Example:  
   ```
      $ CC=cc CXX=CC scons --with-mpi --extrapath=$BOOST_ROOT,$TIFFDIR
   ```

*NOTE: For those users who are member of the 'als' group at NERSC, an installation of `libtiff` is available under `/project/projectdirs/als/local/tiff-4.0.6`.*
4. On successful build, the binary will be generated in the `bin` directory, called `hipgisaxs`.

##### Cori Phase 1 (Cray XC40, Intel Haswell)
1. All required software, except `libtiff` are available as modules on Cori. An example set of modules you can load are given in the `build/modules.cori` file. You could just `source` this file:  
   ```
      $ source build/modules.cori
   ```

2. You will need to install `libtiff`. Please refer to the required software section above.
3. Build command is same as for the Edison system (see above), given in the `build/build-edison.sh`. If needed, make sure the paths are correctly set, including your installation of `libtiff`. Since Cori requires cross compilation for its compute nodes, make sure the `CC` and `CXX` environment variables are also set. Example:

    ```
      $ CC=cc CXX=CC scons --with-mpi --extrapath=$BOOST_ROOT,$TIFFDIR
    ```

*NOTE: For those users who are member of the `als` group at NERSC, an installation of `libtiff` is available under `/project/projectdirs/als/local/tiff-4.0.6`.*
4. On successful build, the binary will be generated in the `bin` directory, called `hipgisaxs`.

#### C. [OLCF](https://www.olcf.ornl.gov) Systems

##### 1. Titan (Cray XK7)
1. All required software as modules on Titan, except `libtiff` which is already installed systemwide. An example set of modules you can load are given in the `build/modules.titan` file. You could just `source` this file:  
   ```
      $ source build/modules.titan
   ```

2. An example build command is given in the `build/build-titan.sh`. If needed, make sure the paths are correctly set. Since Titan requires cross compilation for its compute nodes, make sure the `CC` and `CXX` environment variables are also set. Example:  
   ```
      $ CUDA_TOOLKIT_PATH=$CRAY_CUDATOOLKIT_DIR CXX=CC CC=cc scons --with-mpi --with-cuda --extrapath=$BOOST_DIR
   ```

4. On successful build, the binary will be generated in the `bin` directory, called `hipgisaxs`.

#### D. For users on Andromeda/Bragg servers

1. All the required softwares are already available on these systems. If you want to use your own installation of any of the required softwares, make sure you use its corresponding paths in the following.
2. An example build command is given in the `build/build-bragg.sh` or `build/build-andromeda.sh`. If needed, make sure the paths are correctly set.
3. Once all the required software are available, use the `scons` command to build the binary. Example:  
    ```
      $ scons --extrapath=/usr/local/cuda --with-mpi --with-cuda
    ```

4. On successful build, the binary will be generated in the `bin` directory, called `hipgisaxs`.

### Installing Software Dependencies

#### A. TIFF library

1. The source for the TIFF library can be obtained from http://remotesensing.org/libtiff.
2. Use the `configure` script to generate the required build files. C++ support should be enabled. Additionally, disable JPEG support unless you are willing to install JPEG libraries as well. Example:    
    ```
      $ ./configure --prefix=<my_install_path> --disable-zlib --disable-jpeg --enable-cxx
    ```
3. Compile the source with `make`:  
    ```
      $ make
    ```
4. Install the library at your specified path `<my_install_path>`:  
    ```
      $ make install
    ```
