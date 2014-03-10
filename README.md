# HipGISAXS Version 1.0: QUICK STARTER GUIDE #
Bug (what bug?) reporting: Email: bug-submit@saxs-waxs-gpu.dhcp.lbl.gov and asarje@lbl.gov

## TABLE OF CONTENTS ##
  1. LICENSING
  2. SUPPORTED PLATFORMS
    A. Operating System Platforms
    B. System Hardware (Compute Environment)
    C. System Hardware (Processor Architectures)
  3. SOFTWARE PRE-REQUISITES
  4. HipGISAXS DIRECTORY LAYOUT
  5. TO BUILD THE APPLICATION AND LIBRARY
  6. TO USE THE HIPGISAXS LIBRARY IN YOUR OWN APPLICATION
  7. TO RUN THE APPLICATION
    A. Interactively on Dirac (NERSC)
    B. Submit batch job for multiple GPU nodes on Dirac (NERSC)
    C. Interactively on a generic Linux machine equipped with GPU
  8. INPUTS
  9. APPENDIX
    A. SPECIAL BUILD INSTRUCTIONS FOR SPECIFIC SYSTEMS
      i. On Carver/Dirac at NERSC
      ii. On Hopper/Edison at NERSC
      iii. On Titan at OLCF
      iv. On Stampede at TACC
      v. On Mira at ALCF
      vi. On a generic Linux system
      vii. For certain users, on ALS GPU servers


## LICENSING ##
The HipGISAXS software is only available to be downloaded and used by employees
of academic research institutions, not-for-profit research laboratories, or
governmental research facilities. Please read the accompanying LICENSE file
before downloading the software. By downloading the software, you are agreeing
to be bound by the terms of this Non-Commercial End User License Agreement.
This licensing information is also provided as a PDF file along with the documentation.


## SUPPORTED PLATFORMS
HipGISAXS has been tested on, and supports, the following.

### A. Operating System Platforms
1. GNU/Linux x86\_64: Ubuntu, Red Hat Linux, SUSE Linux, Cray Linux Environment (XT5, XE5, XK6, XE6, XC30).
2. Darwin x86\_64: Mac OS X (Lion, Mountain Lion).
3. You could try HipGISAXS on any UNIX based OS, but it may or may not support it.
   Please let us know about your platform (except Windows!) so that we can include it in our list.

### B. System Hardware (Compute Environment)
1. Clusters/Supercomputers equipped with Nvidia GPUs as accelerators on each node when using the GPU version.
2. Generic desktop equipped with an Nvidia GPU when using the GPU version.
3. Clusters/Supercomputers equipped with Intel Phi coprocessors on each node when using MIC version.
4. Any Intel CPU equipped with Intel Phi coprocessor when using MIC version.
5. Any generic CPU (multi-cores), without accelerators.

### C. System Hardware (Processor Architectures)
1. Intel/AMD processors (including multi-cores).
   We have not extensively tested on various types of processors yet.
2. Fermi or Kepler architecture based Nvidia GPUs.
   Compute capability of these GPUs >= 2.0.
3. Intel MIC architecture based Intel Phi coprocessors.
   Requires intel compilers, version >= 13.1.0.146.


## SOFTWARE PRE-REQUISITES
This software uses several third-party libraries, and they need to be installed and available in order to compile and run HipGISAXS.
For ease, if possible use the installations already available on your system, if any. Alternatively, download and install them yourself.
The following are the dependencies of this software:
1. GNU C/C++ compilers, version >= 4.3 and <= 4.7 (4.6 recommended), for CPU and GPU versions, OR
   Intel C/C++ compilers, version >= 13.1.0.146, necessary for MIC version. Can be used for other versions, but has not been tested. Intel compilers will work best with GNU compatibility 4.6.
2. Nvidia CUDA version > 4.x. (>= 5.0 recommended).
   CUDA can be obtained from: http://developer.nvidia.com/cuda/cuda-downloads
   NOTE: CUDA is NOT required if compiling for CPU-only or MIC versions.
3. GNU (or Intel depending on your compilers) compiled OpenMPI version > 1.4.4.
   OpenMPI can be obtained from: http://www.open-mpi.org/software
   NOTE: Currently, even if you do not plan to use MPI, the compilation process still needs MPI.
   		 This requirement will be removed in future.
   Alternative implementations of MPI can also be used, such as MPICH2 and MVAPICH.
4. Boost and numeric extension to Boost GIL.
   Boost can be obtained from: http://www.boost.org
   Numeric extension to Boost GIL can be obtained from the following (file name 'numeric.zip'):
     http://sourceforge.net/adobe/genimglib/wiki/Downloads
     http://gil-contributions.googlecode.com/svn/trunk
   After downloading, unzip the file, generating a directory named 'numeric'.
   Copy this directory under boost/gil/extensions/ at your Boost installation location.
5. Parallel HDF5 library.
   HDF5 can be obtained from: http://www.hdfgroup.org/downloads
   NOTE: HDF5 depends on zlib and szip.
   zlib can be obtained from: http://www.zlib.net
   szip can be obtained from: http://www.hdfgroup.org/doc\_resource/SZIP
6. Tiff Library (libtiff).
   Tiff library can be obtained from: http://www.libtiff.org
7. Scons version >= 2.0, for installation. It can be obtained from http://www.scons.org


## HipGISAXS DIRECTORY LAYOUT
   hipgisaxs
       |- README     : Duh!
       |- LICENSE    : This contains all the licensing information for HipGISAXS.
       |- SConstruct : Scons file for installation.
       |- SConscript : Scons file for installation.
       |- bin        : This contains HipGISAXS binaries generated by compilation.
       |- build      : A few Makefiles for various systems are provided in this.
       |- data       : This provides some sample input shape definition files in HDF5 format.
       |- doc        : This contains detailed documentation of HipGISAXS.
       |- inputs     : This gives some sample input files to HipGISAXS program in HiG format.
       |- lib        : The HipGISAXS library, libhipgisaxs.a, is generated here.
       |- man        : This has the man pages.
       |- obj        : All the object files generated during build are stored here.
       |- samples    : This contains some compilable samples on how to use the HipGISAXS library.
       |- src        : This is the main source directory containing the whole code.
                         It contains many subdirectories.


## TO BUILD THE APPLICATION AND LIBRARY
To build the HipGISAXS application binary and static library, use "scons". Make sure you are passing paths to all the dependencies through "--extrapath=" option of scons:
  $ scons --extrapath=<path1>,<path2>,<etc>
The generated binary, "hipgisaxs", will be located in the "bin" directory.
To enable GPU support, use "--with-cuda" option:
  $ scons --extrapath=<path1>,<path2>,<etc> --with-gpu
To enable MPI support, use "--with-mpi" option:
  $ scons --extrapath=<path1>,<path2>,<etc> --with-gpu --with-mpi
or,
  $ scons --extrapath=<path1>,<path2>,<etc> --with-mpi
The generated binary, "hipgisaxs", will be located in the "bin" directory.
The generated static library, "libhipgisaxs.a", will be located in the "lib" directory.

... and you are done. Go enjoy HipGISAXS!

NOTE: See Appendix at the end of this file for more detailed and customized building information.


## TO USE THE HIPGISAXS LIBRARY IN YOUR OWN APPLICATION
   Please refer to the example provided in the "samples" directory.
   It contains a simple code which uses the hipgisaxs library.
   Its corresponsing Makefile is also provided.


## TO RUN THE APPLICATION
   
### A. Interactively on Dirac (NERSC)
1. Unload the default PGI modules:
    $ module unload pgi openmpi
2. Load the required modules:
    $ module load openmpi-gnu gcc/4.5.2 cuda/4.2
    $ module load szip zlib
    $ module load hdf5-parallel/1.8.3-gnu
3. Request a GPU node:
    $ qsub -I -V -q dirac_int -l nodes=1:ppn=8:fermi
4. Make sure loaded modules are correct.
   By default, PGI version of openMPI might be loaded, and needs to be unloaded:
    $ module unload openmpi
5. Move to code directory:
    $ cd $PBS_O_WORKDIR
6. Execute on a single node:
    Usage: ./bin/hipgisaxs <input-file-in-HiG>
   Example:
    $ ./bin/hipgisaxs inputs/test.27.hig
7. For more details on running interactive jobs on Dirac, please refer to:
    http://www.nersc.gov/users/computational-systems/dirac/running-jobs/interactive

### B. Submit batch job for multiple GPU nodes on Dirac (NERSC)
1. Create a job script. Example:
     $ cat script.pbs
     #PBS -q dirac_special
     #PBS -l nodes=12:ppn=1:fermi
     #PBS -l walltime=03:00:00
     #PBS -A gpgpu
     #PBS -N opv_new1.12
     #PBS -e opv_new1.12.$PBS_JOBID.err
     #PBS -o opv_new1.12.$PBS_JOBID.out
     #PBS -V
     cd $PBS_O_WORKDIR
     module unload pgi openmpi
     module load openmpi-gnu/1.4.5 gcc/4.5.4 cuda/4.2
     module load szip zlib
     module load hdf5-parallel/1.8.3-gnu
     export PATH=/global/homes/a/asarje/local/tiff-4.0.2/bin:$PATH
     export LD_LIBRARY_PATH=/global/homes/a/asarje/local/tiff-4.0.2/lib:$LD_LIBRARY_PATH
     mpirun -np 4 ./bin/hipgisaxs inputs/test.27.hig
2. In the submission script:
   Make sure the modules are loaded correctly (see example above.)
   For the execution command, use:
     mpitun -np <nodes> ./bin/hipgisaxs <input-file-in-HiG>
   where, 'nodes' is number of GPU nodes.
3. Submit the script:
    $ qsub script.pbs
4. For detailed information on writing and submitting job scripts for Dirac, please refer to:
    http://www.nersc.gov/users/computational-systems/dirac/running-jobs/batch

### C. Interactively on a generic Linux machine equipped with GPU (including saxs-waxs-gpu)
1. Make sure all the paths (data, output) and other variables are correctly set in the input HiG file.
2. Execute the binary:
    Usage: ./bin/hipgisaxs <input-file-in-HiG>
   Example:
    $ ./bin/hipgisaxs inputs/test.27.hig


## INPUTS
The HipGISAXS binary takes as input a file in HiG format.
Please refer to detailed HipGISAXS documentation for details of the HiG format.
A few sample input files are located in the directory "inputs", with extensions ".hig".
Update the input file as needed.
The main components in the input file to update are the following:
1. Shape name defines the input filename containing triangulated shape data.
   It should point to the correct location of the file (relative path):
       ...
       shape = {
    	   ...
    	   name = "data/Shape_27.hd5"
       } ...
   This file needs to be in HDF5 format.
   Some sample shape files are provided in the directory "data", with extensions ".hd5".
2. Output locatation needs to be defined as pathprefix and runname in the computation object:
      ...
      computation = {
    	  ...
    	  pathprefix = ".",
    	  runname = "27",
    	  ...
      } ...
   The pathprefix is a relative path to a directory.
   The runname is appended with a timestamp and a directory by this resulting name
    is created in the directory specified by pathprefix.
   All generated output files are stored in this directory.
3. Resolution alters the final image resolution, and also affects the run time:
      ...
      computation = {
    	  ...
    	  resolution = [ 0.5 0.5 ],
    	  ...
      } ...
   Lower resolutions execute faster, and higher slower.
   Modify it as needed.



# APPENDIX

## SPECIAL BUILD INSTRUCTIONS FOR SPECIFIC SYSTEMS (STALE INFO, TO BE UPDATED SOON)

### A. On Carver/Dirac at NERSC
1. Some makefiles for different systems are included in the directory "build".
   Replace the current "Makefile" with the one for Carver/Dirac (and name it as "Makefile"):
    $ cp build/Makefile.dirac Makefile
2. Unload the default PGI modules:
    $ module unload pgi openmpi
3. Load the required modules:
    $ module load openmpi-gnu gcc/4.5.2 cuda
    $ module load szip zlib
    $ module load hdf5-parallel/1.8.3-gnu
4. The Boost module available on Carver does NOT contain the GIL numeric extension.
   Please install and use your own copy.
5. Edit Makefile to specify the correct paths ("base directories" section) to the various libraries.
   Example:
    $ cat Makefile
    ...
    ## base directories
    BOOST_DIR = /global/homes/a/asarje/local/boost_1_49_0
    MPI_DIR = 
    CUDA_DIR = /usr/common/usg/cuda/5.0
    HDF5_DIR = /global/homes/a/asarje/local/hdf5-1.8.8-gnu/parallel
    TIFF_DIR = /global/homes/a/asarje/local/tiff-4.0.2
    Z_DIR = $(ZLIB_DIR)     # an environment variable set by loading zlib module
    SZ_DIR = $(SZIP_DIR)    # an environment variable set by loading szip module
    ...
6. Build the code from within the code directory:
    $ make clean
    $ make [or, make library]
   This will generate the binary in directory "bin".
   When building the library, it will be generated into the "lib" directory.
   All the intermediate generated object files are in the directory "obj". They can be removed if wanted.

### B. On Hopper/Edison at NERSC
   Instructions coming soon. Meanwhile go try yourself!

### C. On Titan at OLCF
   Instructions coming soon. Meanwhile go try yourself!

### D. On Stampede at TACC
   Instructions coming soon. Meanwhile go try yourself!

### E. On Mira at ALCF
   Instructions coming soon. Meanwhile go try yourself!

### F. On a generic Linux system
1. Make sure all above prerequisites are available.
2. Make sure all system environment variables are set accordingly to include the prerequisites.
3. Set all paths correctly in the sample Makefile, and edit as needed (as above).
4. Build the library from within the code directory:
    $ make clean
    $ make library

### G. For certain users, on saxs-waxs-gpu/saxs-waxs-gpu2/andromeda servers
1. Replace the current "Makefile" with the one for saxs-waxs-gpu (and name it as "Makefile"):
    $ cp build/Makefile.saxs1 Makefile		OR
    $ cp build/Makefile.saxs2 Makefile		OR
    $ cp build/Makefile.andromeda Makefile
2. All "base directories" in the provided Makefile are set to the correct locations.
   If you want to use your own installation of any of the required softwares, please edit its corresponding entry in the Makefile.
   Example:
    $ cat Makefile
    ...
    ## base directories
    BOOST_DIR = /usr/local/boost_1_45_0
    MPI_DIR = /usr/local
    CUDA_DIR = /usr/local/cuda
    HDF5_DIR = /home/asarje/local/hdf5-1.8.8-gnu/parallel
    Z_DIR = /root/zlib-1.2.7
    SZ_DIR = /root/szip-2.1
    TIFF_LIB_DIR = /usr/local
    ...
3. Build the code from within the code's main directory:
    $ make clean
    $ make [or, make library]
   This will generate the binary in directory "bin".
   When building the library, it will be generated into the "lib" directory.
   All the generated object files are in the directory "obj".