import os
import sys
import textwrap
import copy
import re

EnsureSConsVersion(1, 1, 0)

## detect which platform

platform = os.sys.platform
nix = False
linux = False
darwin = False
windows = False
freebsd = False
openbsd = False
solaris = False

if "darwin" == platform:
    darwin = True
    platform = "osx"
elif platform.startswith("linux"):
    linux = True
    platform = "linux"
elif "sunos5" == platform:
    solaris = True
elif platform.startswith("freebsd"):
    freebsd = True
elif platform.startswith("openbsd"):
    openbsd = True
elif "win32" == platform:
    windows = True
else:
    print( "No special config for [" + platform + "] which probably means it won't work" )

nix = not windows

if windows:
    print("There is no Windows version of HipGISAXS :P")
    Exit(0)


## some global variables

DEFAULT_INSTALL_DIR = "."
toolchain_gcc = "GNU"
toolchain_intel = "Intel"
toolchain_pgi = "PGI"


## function definitions

def add_option(name, help, nargs, contributesToVariantDir, dest = None, default = None,
                type = "string", choices = None, metavar = None):
    if dest is None:
        dest = name
    if type == 'choice' and not metavar:
        metavar = '[' + '|'.join(choices) + ']'
    AddOption("--" + name, dest = dest, type = type, nargs = nargs, action = "store",
                choices = choices, default = default, metavar = metavar, help = help)
    options[name] = {   "help"   : help,
                        "nargs"  : nargs ,
                        "dest"   : dest,
                        "default": default
                    }

def get_option(name):
    return GetOption(name)

def _has_option(name):
    x = get_option(name)
    if x is None:
        return False
    if x == False:
        return False
    if x == "":
        return False
    return True

def printLocalInfo():
    import sys, SCons
    print("scons version: " + SCons.__version__)
    print("python version: " + " ".join([`i` for i in sys.version_info]))

def add_paths(s):
    for x in s.split(","):
        env.Append(CPPPATH = [x + "/include"])
        env.Append(LIBPATH = [x + "/lib64"])
        env.Append(LIBPATH = [x + "/lib"])
        path = env['ENV']['PATH'] + ":" + x + "/bin"
        env['ENV']['PATH'] = path

def add_include_path(p):
  env.Append(CPPPATH = [ p ])

def add_library_path(p):
  env.Append(LIBPATH = [ p ])

def add_libs(s):
    for x in s.split(","):
        env.Append(LIBS = [x])

def add_ld_library_path(env, s):
    env['ENV']['LD_LIBRARY_PATH'] = s
    for x in s.split(":"):
        env.Append(LIBPATH = [x])

## setup configuration tests which can be done
def setup_configuration_tests(conf):
    def FindSysLibDep(context, name, libs, **kwargs):
        var = "LIBDEPS_" + name.upper() + "_SYSLIBDEP"
        var = var.replace("++", "PP")
        kwargs['autoadd'] = False
        for lib in libs:
            result = context.sconf.CheckLib(lib, **kwargs)
            context.did_show_result = 1
            if result:
                context.env[var] = lib
                return context.Result(result)
        context.env[var] = name
        return context.Result(result)

    conf.AddTest('FindSysLibDep', FindSysLibDep)

    ## identify the toolchain
    def check_for_toolchain(context, toolchain, lang_name, compiler_var, source_suffix):
        test_bodies = {
                toolchain_gcc: ("""
                                #if !defined(__GNUC__) && !defined(__GNUG__)
                                #error
                                #endif
                                """),
                toolchain_intel: ("""
                                #if !defined(__INTEL_COMPILER) && !defined(__ICC)
                                #error
                                #endif
                                """),
                toolchain_pgi: ("""
                                #if !defined(__PGI)
                                #error
                                #endif
                                """)
                }
        print_tuple = (lang_name, context.env[compiler_var], toolchain)
        context.Message('Checking if %s compiler "%s" is %s ... ' % print_tuple)
        # Strip indentation from the test body to ensure that the newline at the end of the
        # endif is the last character in the file (rather than a line of spaces with no
        # newline), and that all of the preprocessor directives start at column zero. Both of
        # these issues can trip up older toolchains.
        test_body = textwrap.dedent(test_bodies[toolchain])
        result = context.TryCompile(test_body, source_suffix)
        context.Result(result)
        return result

    conf.AddTest('check_for_toolchain', check_for_toolchain)

    #def check_for_accelerator(context, toolchain, lang_name, compiler_var, source_suffix):
    #   test_bodies = {
    #           accelerator_: ("""
    #                           #if !defined(__GNUC__) && !defined(__GNUG__)
    #                           #error
    #                           #endif
    #                           """),
    #           toolchain_intel: ("""
    #                           #if !defined(__INTEL_COMPILER) && !defined(__ICC)
    #                           #error
    #                           #endif
    #                           """),
    #           toolchain_pgi: ("""
    #                           #if !defined(__PGI)
    #                           #error
    #                           #endif
    #                           """)
    #           }
    #   print_tuple = (lang_name, context.env[compiler_var], toolchain)
    #   print(print_tuple)
    #   context.Message('Checking if CUDA works ... ')
    #   test_body = textwrap.dedent(test_bodies[toolchain])
    #   result = context.TryCompile(test_body, source_suffix)
    #   context.Result(result)
    #   return result

#if defined(__clang__)  /* Clang/LLVM. */
#elif defined(__ICC) || defined(__INTEL_COMPILER)   /* Intel ICC/ICPC. */
#elif defined(__GNUC__) || defined(__GNUG__)    /* GNU GCC/G++. */
#elif defined(__HP_cc) || defined(__HP_aCC) /* Hewlett-Packard C/aC++. */
#elif defined(__IBMC__) || defined(__IBMCPP__)  /* IBM XL C/C++. */
#elif defined(_MSC_VER) /* Microsoft Visual Studio. */
#elif defined(__PGI)    /* Portland Group PGCC/PGCPP. */
#elif defined(__SUNPRO_C) || defined(__SUNPRO_CC)   /* Oracle Solaris Studio. */
#endif

def do_configure(myenv):
    def add_flag_if_supported(env, tool, extension, flag, **mutation):
        def check_flag_test(context, tool, extension, flag):
            test_body = ""
            context.Message("Checking if %s compiler supports %s ... " % (tool, flag))
            ret = context.TryCompile(test_body, extension)
            context.Result(ret)
            return ret

        test_mutation = mutation
        if using_gcc():
            test_mutation = copy.deepcopy(mutation)
            for key in test_mutation.keys():
                test_flags = test_mutation[key]
                for test_flag in test_flags:
                    if test_flag.startswith("-Wno-") and not test_flag.startswith("-Wno-error="):
                        test_flags.append(re.sub("^-Wno-", "-W", test_flag))

        cloned = env.Clone()
        cloned.Append(**test_mutation)

        cloned.Append(CCFLAGS = ['-Werror'])
        conf = Configure(cloned, help = False, clean = False, custom_tests = {
                        'CheckFlag': lambda(ctx): check_flag_test(ctx, tool, extension, flag)
            })
        available = conf.CheckFlag()
        conf.Finish()
        if available:
            env.Append(**mutation)
        return available

    def add_to_cflags_if_supported(env, flag):
        return add_flag_if_supported(env, 'C', '.c', flag, CFLAGS = [flag])

    def add_to_ccflags_if_supported(env, flag):
        return add_flag_if_supported(env, 'C', '.c', flag, CCFLAGS = [flag])
    
    def add_to_cxxflags_if_supported(env, flag):
        return add_flag_if_supported(env, 'C++', '.cpp', flag, CXXFLAGS = [flag])
    
    def add_to_linkflags_if_supported(env, flag):
        return add_flag_if_supported(env, 'C', '.c', flag, LINKFLAGS = [flag])
    
    ## configure

    conf = Configure(myenv, help = False, clean = False, custom_tests = {
    #                                       'check_for_toolchain': check_for_toolchain,
    #                                       'FindSysLibDep': FindSysLibDep
                                            })
    #print dir(conf)

    ## setup configuration tests
    setup_configuration_tests(conf)
    
    if 'CXX' in os.environ:
        conf.env.Replace(CXX = os.environ['CXX'])
        print(">> Using C++ compiler " + os.environ['CXX'])
    if 'CC' in os.environ:
        conf.env.Replace(CC = os.environ['CC'])
        print(">> Using C compiler " + os.environ['CC'])

    ## select compiler toolchain
    toolchain = None
    have_toolchain = lambda: toolchain != None
    using_gcc = lambda: toolchain == toolchain_gcc
    using_intel = lambda: toolchain == toolchain_intel
    using_pgi = lambda: toolchain == toolchain_pgi
    toolchain_search_sequence = [toolchain_intel, toolchain_pgi, toolchain_gcc] ## note GNU is last
    for candidate in toolchain_search_sequence:
        if conf.check_for_toolchain(candidate, "C++", "CXX", ".cpp"):
            toolchain = candidate
            break
    #if not get_option('clean'):
    if not have_toolchain():
        print("error: could not identify toolchain on your system")
        Exit(1)
    if not conf.check_for_toolchain(toolchain, "C", "CC", ".c"):
        print("error: C toolchain does not match the identified C++ toolchain")
        Exit(1)
    conf.env['TOOLCHAIN'] = toolchain

    if not conf.CheckCHeader('assert.h'):
        print("error: C header 'assert.h' missing in your compiler installation!")
        Exit(1)

    ## check if compilers work
    if not conf.CheckCXX():
        print("error: your C++ compiler '%s' doesn't seem to work!" % conf.env["CXX"])
        Exit(1)
    if not conf.CheckCC():
        print("error: your C compiler '%s' doesn't seem to work!" % conf.env["CC"])
        Exit(1)

    ## check for functions
    if not conf.CheckFunc('printf'):
        print('error: your compiler and/or environment is not correctly set.')
        Exit(0)

    ## check for headers
    if not conf.CheckCXXHeader('iostream'):
        print "error: you need 'iostream' to build"
        Exit(1)
    if not conf.CheckCXXHeader('vector'):
        print("error: you need 'std::vector' to build")
        Exit(1)
    if not conf.CheckCXXHeader("boost/filesystem/operations.hpp"):
        print( "can't find boost headers" )
        Exit(1)
    if using_parallel_hdf5:
        if not conf.CheckCHeader('hdf5.h'):
            print("error: you need HDF5 header to build")
            Exit(1)
    conf.env.Append(CPPDEFINES = [("BOOST_THREAD_VERSION", "2")])
    #for b in boost_libs:
    #   l = "boost_" + b
    #   conf.FindSysLibDep(l, [l + boost_c + "-mt" + boost_v,
    #                           l + boost_c + boost_v ], language='C++')

#   if not conf.CheckLib('cuda') or not conf.CheckLib('cudart'):
#       print("warning: CUDA not found. proceeding without GPU support.")
#   else:
#       print("CUDA found. Proceeding with GPU support.")

    if not using_tau:
      ## check for types
      if conf.CheckType('float_t') or conf.CheckType('complex_t'):
          print("warning: Some types used are already defined. Overriding them.")

    myenv = conf.Finish()

    ## add flags
    if using_gcc() or using_pgi():
        add_to_ccflags_if_supported(myenv, '-Wno-unused-local-typedefs')
    if not add_to_cxxflags_if_supported(myenv, '-std=c++11'):
        if not add_to_cxxflags_if_supported(myenv, '-std=c++0x'):
            print("C++11 or C++0x mode is necessary. Cannot find flag to enable it.")
            Exit(1)
    if not add_to_cflags_if_supported(myenv, '-std=c99'):
        print('C++11 enabled for C++ sources, but failed to enable C99 for C sources')
    if using_gcc():
        add_to_ccflags_if_supported(myenv, '-fopenmp')
    if using_intel():
        add_to_ccflags_if_supported(myenv, '-qopenmp')
        add_to_linkflags_if_supported(myenv, '-qopenmp')

    return myenv


def do_post_configure(myenv):
    conf = Configure(myenv, help = False, clean = False)

    ## setup configuration tests
    setup_configuration_tests(conf)
    
    ## check for libraries
    for libs in all_libs:
        conf.FindSysLibDep(libs, [libs], language='C++')

    myenv = conf.Finish()
    return myenv



################
## main stuff ##
################

## options and options

options = { }

add_option("mute" , "reduce screen noise", 0, False, default = 0)
add_option("prefix", "installation prefix", 1, False, default = DEFAULT_INSTALL_DIR)
# compiler stuff
add_option("cxx", "compiler to use for c++", 1, True)
add_option("cc", "compiler to use for c", 1, True)
add_option("ld", "linker to use", 1, True)
add_option("cpppath", "Include path if you have headers in a nonstandard directory", 1, True)
add_option("libpath", "Library path if you have libraries in a nonstandard directory", 1, True)
add_option("extrapath", "comma separated list of additional paths", 1, True)
add_option("extralib", "comma separated list of additional libraries", 1, True)
## debug options
add_option("release", "release build, full optimizations, no debug", 0, True, default = 0)
add_option("dbg", "debug build, no optimization", 0, True, "debug_build", default = 0)
add_option("dbgo", "optimized build with debug symbols", 0, True, "opt_debug_build", default = 0)
add_option("dbgl", "debug build, no optimization, debug logging", 0, True,
            "debug_build_and_log", default = 0)
## optional packages and support
add_option("with-cuda", "Enable GPU support", 0, False)
add_option("with-mic", "Enable Intel MIC support", 0, False)
add_option("with-mpi", "Enable MPI parallelization", 0, False)
add_option("with-parallel-hdf5", "Enable parallel HDF5 data format support", 0, False)
add_option("with-papi", "Enable PAPI profiling", 0, False)
add_option("use-single", "Use single precision floating point arithmatic", 0, False)
## other stuff
add_option("detail-time", "Output detailed timings", 0, False)
add_option("detail-mem", "Output detailed memory usage", 0, False)
add_option("verbose", "Be verbose", 0, False)
add_option("with-tau", "Compile with TAU for profiling", 0, False)
add_option("with-yaml", "Use YAML input files.", 0, False)

printLocalInfo()

# do not build if help is called
if GetOption('help'):
    Return()

linux64 = True
print("Building for 64-bit")
msarch = "amd64"

## get the optional command line parameters
release_build = _has_option("release")
using_debug = _has_option("debug_build") or _has_option("debug_build_and_log")
using_opt_debug = _has_option("opt_debug_build")
if release_build and using_debug:
    print("error: release build cannot have debug enabled. make a choice first.")
    Exit(1)
using_cuda = _has_option("with-cuda")
using_mic = _has_option("with-mic")      ## disabling MIC version
using_mpi = _has_option("with-mpi")
using_parallel_hdf5 = _has_option("with-parallel-hdf5")
using_papi = _has_option("with-papi")
using_single = _has_option("use-single")
using_tau = _has_option("with-tau")
using_yaml = _has_option("with-yaml")

if not using_mpi and using_parallel_hdf5:
    print("error: to enable parallel HDF5 support, you need to enable MPI as well.")
    Exit(1)

if using_cuda and using_mic:
    print("error: currently GPU and MIC are not supported simultaneously. select one of them.")
    Exit(1)

#variant_dir = ".build"
variant_dir = "obj"

#########################
## set the environment ##
#########################

if using_mpi:
    ## general
    CXXCOMPILER = "mpicxx"
    CCCOMPILER = "mpicc"
else:
    CXXCOMPILER = "c++"
    CCCOMPILER = "cc"
if using_tau:
  CXXCOMPILER = "tau_cxx.sh"
  CCCOMPILER = "tau_cc.sh"

if 'CXX' in os.environ: CXXCOMPILER = os.environ['CXX']
if 'CC' in os.environ: CCCOMPILER = os.environ['CC']

# initialize the environment with gathered stuff so far
env = Environment(BUILD_DIR=variant_dir,
                    ENV = os.environ,
                    CC = CCCOMPILER,
                    CXX = CXXCOMPILER,
                    DIST_ARCHIVE_SUFFIX = '.tar',
                    TARGET_ARCH = msarch,
                    PYSYSPLATFORM = os.sys.platform)
SCONSIGN_FILE = ".scons-signatures"
env.SConsignFile(SCONSIGN_FILE)

## the decider
env.Decider('MD5-timestamp')

## if muted
if _has_option("mute"):
    env.Append(CCCOMSTR = "Compiling $TARGET")
    env.Append(CXXCOMSTR = env["CCCOMSTR"])
    env.Append(SHCCCOMSTR = env["CCCOMSTR"])
    env.Append(SHCXXCOMSTR = env["SHCCCOMSTR"])
    env.Append(LINKCOMSTR = "Linking $TARGET")
    env.Append(SHLINKCOMSTR = env["LINKCOMSTR"])
    env.Append(ARCOMSTR = "Generating library $TARGET")

if "uname" in dir(os):
    processor = os.uname()[4]
else:
    processor = "x86_64"
env['PROCESSOR_ARCHITECTURE'] = processor
print("Detected processor architecture: %s" % processor)
env.Append(CCFLAGS = "-m64")

install_dir = DEFAULT_INSTALL_DIR
if _has_option("prefix"):
    install_dir = get_option("prefix")
env['INSTALL_DIR'] = install_dir


if _has_option("cpppath"):
    env.Append(CPPPATH =  get_option("cpppath"))

## extra paths and libs
working_dir = env['ENV']['PWD']
env.Append(CPPPATH = [working_dir + "/include"])
if _has_option("extrapath"):
    add_paths(get_option("extrapath"))
if _has_option("extralib"):
    add_libs(get_option("extralib"))


nix_lib_prefix = "lib"
if processor == "x86_64":
    linux64 = True
    nix_lib_prefix = "lib64"
#    env.Append(LIBPATH = ["/usr/lib64",
#                            "/usr/lib",
#                            "/lib64",
#                            "/lib",
#                            "/usr/local/lib64",
#                            "/usr/local/lib"])
    #env.Append(CPPPATH = [ '/usr/include/x86_64-linux-gnu/c++/4.8' ])
env['NIX_LIB_DIR'] = nix_lib_prefix

#add_ld_library_path(env, os.environ['LD_LIBRARY_PATH'])

if _has_option("cpppath"):
    env.Append(CPPPATH = [get_option("cpppath")])

#print env['LIBPATH']
#print env['ENV']

use_mkl = False

## call configure
if not get_option('clean'):
    env = do_configure(env)

    ## MKL with intel only (Edison)
    if env['TOOLCHAIN'] == toolchain_intel:
      use_mkl = True
      if 'MKL_INC' in env['ENV']:
        add_include_path(env['ENV']['MKL_INC'])
      if 'MKL_LIBDIR' in env['ENV']:
        add_library_path(env['ENV']['MKL_LIBDIR'])

    ## stuff for TAU
    if using_tau:
      if 'TAU_MAKEFILE' not in env['ENV']:
        print "error: please specify TAU_MAKEFILE environment variable to enable use of TAU"
        sys.exit(1)

    using_accelerator = None
    ## required libs
    boost_libs = ["boost_system", "boost_filesystem", "boost_timer", "boost_chrono"]
    parallel_hdf5_libs = ["hdf5", "z"]
    tiff_libs = ["tiff"]
    mkl_libs = []
    if use_mkl:
      mkl_libs = [] #["mkl_intel_lp64", "mkl_sequential", "mkl_core"]
    if env['TOOLCHAIN'] == toolchain_intel:
      other_libs = ["m", "stdc++"]
      if using_mic:
        other_libs.append("pfm")
    elif platform == "osx":
        other_libs = [ "m" ]
    else:
      other_libs = ["m", "gomp"]
    ## optional libs
    mpi_libs = []
    #mpi_libs = ["mpi_cxx", "mpi"]
    #gpu_libs = ["cudart", "cufft", "cudadevrt", "nvToolsExt"]
    gpu_libs = ["cudart"]
    papi_libs = ["papi"]
    ## required flags
    detail_flags = []
    if _has_option("detail-time"): detail_flags += ['TIME_DETAIL_1', 'TIME_DETAIL_2']
    if _has_option("detail-mem"): detail_flags += ['MEM_DETAIL']
    if _has_option("verbose"): detail_flags += ['FF_VERBOSE', 'SF_VERBOSE']
    ## optional flags
    gpu_flags = ['USE_GPU', 'GPUR', 'KERNEL2', 'FF_ANA_GPU', 'FF_NUM_GPU', 'FF_NUM_GPU_FUSED', 'SF_GPU']
    mic_flags = [] #['USE_MIC', 'FF_MIC_OPT', 'FF_NUM_MIC_SWAP', 'FF_NUM_MIC_KB']
    cpu_flags = ['FF_NUM_CPU_FUSED', 'FF_CPU_OPT']
    ## choose at most one of the following
    #if use_mkl: cpu_flags += [ 'FF_CPU_OPT_MKL' ]
    cpu_flags += [ 'INTEL_AVX', 'FF_CPU_OPT_AVX' ]

    mpi_flags = ['USE_MPI']
    parallel_hdf5_flags = ['USE_PARALLEL_HDF5']
    papi_flags = ['PROFILE_PAPI']

    ioflags = ['FILEIO']

    all_flags = detail_flags + ioflags
    all_libs = boost_libs + tiff_libs + mkl_libs + other_libs

    if using_cuda:
        all_libs += gpu_libs
        all_flags += gpu_flags
        using_accelerator = 'gpu'
    elif using_mic:
        all_flags += mic_flags + cpu_flags
        using_accelerator = 'mic'
        env.Append(CCFLAGS = ["-mmic", "-mt_mpi"])
        env.Append(LINKFLAGS = ["-mmic", "-mt_mpi", "-mkl"])
    else:
        all_flags += cpu_flags

    if using_mpi:
        all_libs += mpi_libs
        all_flags += mpi_flags

    if using_parallel_hdf5:
        all_libs += parallel_hdf5_libs
        all_flags += parallel_hdf5_flags
        env['USE_PARALLEL_HDF5'] = True
    else:
        env['USE_PARALLEL_HDF5'] = False

    if using_papi:
        all_libs += papi_libs
        all_flags += papi_flags

    env.Append(LIBS = all_libs)
    env.Append(CPPDEFINES = all_flags)
    env['ACCELERATOR_TYPE'] = using_accelerator
    if using_accelerator != None:
        print("Enabling use of accelerator: %s" % using_accelerator)

    if not using_mic:
      env.Append(CCFLAGS = ["-march=core-avx2"]) #, "-msse4.1", "-msse4.2", "-mssse3"])
      #env.Append(CCFLAGS = ["-no-vec"])

    if using_debug:
        env.Append(CCFLAGS = ['-g'])
        env.Append(CPPDEFINES = ['DEBUG'])
    elif using_opt_debug:
        odbg_flags = ['-O3', '-g']
        if env['TOOLCHAIN'] == toolchain_intel: odbg_flags += ['-dynamic']
        env.Append(CCFLAGS = odbg_flags)
        env.Append(LINKFLAGS = odbg_flags)
        env.Append(CPPDEFINES = ['NDEBUG'])
    else:
        env.Append(CCFLAGS = ['-O3'])
        env.Append(CPPDEFINES = ['NDEBUG'])

    # use double-precision as defualt
    if not using_single:
        env.Append(CPPDEFINES = ['DOUBLEP'])

    # yaml
    if using_yaml:
      yaml_headers = ['yaml-cpp/yaml.h']
      yaml_libs = ['yaml-cpp']
      env.Append(LIBS = yaml_libs)
      env.Append(CPPDEFINES = ['YAML'])

    ## print stuff
    #for item in sorted(env.Dictionary().items()):
    #   print "++ '%s', value = '%s'" % item
    print("Platform: %s" % sys.platform)

## call post configure
if not get_option('clean'):
    env = do_post_configure(env)

    gpuenv = env.Clone(CC = 'nvcc')
    if using_cuda:
        gpuenv.Tool('cuda')
        ## remove openmp flag from CCFLAGS and LINKFLAGS
        old_ccflags = gpuenv['CCFLAGS']
        ccflags = []
        flags_to_remove = ['-fopenmp', '-openmp', '-Wno-unused-local-typedefs', "-mavx", "-mtune=core-avx2"]
        for flag in old_ccflags:
            if flag not in flags_to_remove:
                ccflags += [ flag ]
        gpuenv.Replace(CCFLAGS = ccflags)
        old_linkflags = gpuenv['LINKFLAGS']
        linkflags = []
        for flag in old_linkflags:
            if flag not in flags_to_remove:
                linkflags += [ flag ]
        gpuenv.Replace(LINKFLAGS = linkflags)
        ## add other flags
        gpuenv.Append(LINKFLAGS = ['-Xlinker', '-Wl,-rpath', '-Xlinker', '-Wl,$CUDA_TOOLKIT_PATH/lib64'])
        gpuenv.Append(LINKFLAGS = ['-Xlinker', '-lgomp'])
        gpuenv.Append(LINKFLAGS = ['-arch=sm_35'])
        gpuenv.Append(LINKFLAGS = ['-dlink'])
    Export('gpuenv')
    Export('env')

    objs, nvobjs, mainobjs = env.SConscript('src/SConscript', duplicate = False)

    if using_cuda:
        nvlibobj = gpuenv.Program('nv_hipgisaxs.o', nvobjs)
        objs += nvlibobj + nvobjs

    env.Library('hipgisaxs', objs)
    env.Program('hipgisaxs', objs + mainobjs)
    env.Install('#/bin', source = 'hipgisaxs')
    env.Install('#/lib', source = 'libhipgisaxs.a')

## to clean everything
Clean('.', '#/' + env['BUILD_DIR'])
Clean('.', '#/bin/hipgisaxs')
Clean('.', '#/lib/libhipgisaxs.a')
Clean('.', env['CONFIGUREDIR'])
Clean('.', env['CONFIGURELOG'])
Clean('.', '#/' + SCONSIGN_FILE + ".dblite")
