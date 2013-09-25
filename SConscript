import os
import sys

EnsureSConsVersion(1, 1, 0)

######## platform - detect which platform
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
elif platform.startswith( "freebsd" ):
	freebsd = True
elif platform.startswith( "openbsd" ):
	openbsd = True
elif "win32" == platform:
	windows = True
else:
	print( "No special config for [" + platform + "] which probably means it won't work" )

nix = not windows

if windows:
	print("There is no Windows version of HipGISAXS :P")
	Exit(0)


DEFAULT_INSTALL_DIR = "."

def add_option(name, help, nargs, contributesToVariantDir, dest = None, default = None,
				type = "string", choices = None, metavar = None):
	if dest is None:
		dest = name
	if type == 'choice' and not metavar:
		metavar = '[' + '|'.join(choices) + ']'
	AddOption("--" + name, dest = dest, type = type, nargs = nargs, action = "store",
				choices = choices, default = default, metavar = metavar, help = help)
	options[name] = {	"help"   : help,
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

def add_extra_paths(s):
	for x in s.split(","):
		env.Append(EXTRACPPPATH = [x + "/include"])
		env.Append(EXTRALIBPATH = [x + "/lib64"])
		env.Append(EXTRALIBPATH = [x + "/lib"])

def add_extra_libs(s):
	for x in s.split(","):
		env.Append(LIBS = [x])

def do_configure(myenv):
	## check if compilers work
	conf = Configure(myenv, help = False)
	if 'CheckCXX' in dir(conf):
		if not conf.CheckCXX():
			print("C++ compiler %s does not work" % conf.env["CXX"])
			Exit(1)
	if myenv["CC"] != myenv["CXX"]) and 'CheckCC' in dir(conf):
		if not conf.CheckCC():
			print("C compiler %s does not work" % conf.env["CC"])
			Exit(1)
	myenv = conf.Finish()

	## identify the toolchain
	toolchain_gcc = "GCC"
	toolchain_intel = "Intel"

	#if defined(__clang__)	/* Clang/LLVM. */
	#elif defined(__ICC) || defined(__INTEL_COMPILER)	/* Intel ICC/ICPC. */
	#elif defined(__GNUC__) || defined(__GNUG__)	/* GNU GCC/G++. */
	#elif defined(__HP_cc) || defined(__HP_aCC)	/* Hewlett-Packard C/aC++. */
	#elif defined(__IBMC__) || defined(__IBMCPP__)	/* IBM XL C/C++. */
	#elif defined(_MSC_VER)	/* Microsoft Visual Studio. */
	#elif defined(__PGI)	/* Portland Group PGCC/PGCPP. */
	#elif defined(__SUNPRO_C) || defined(__SUNPRO_CC)	/* Oracle Solaris Studio. */
	#endif
	def check_for_toolchain(context, toolchain, lang_name, compiler_var, source_suffix):
		test_bodies = {
				toolchain_gcc: ("""
								#if !defined(__GNUC__) || defined(__clang__)
								#error
								#endif
								"""),
				toolchain_intel: ("""
								#if !defined(__INTEL_COMPILER) || !defined(__ICC)
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

	conf = Configure(myenv, help = False, custom_tests = {'CheckForToolChain': check_for_toolchain})
	toolchain = None
	have_toolchain = lambda: toolchain != None
	using_gcc = lambda: toolchain == toolchain_gcc
	using_intel = lambda: toolchain == toolchain_intel
	using_pgi = lambda: toolchain == toolchain_pgi

	toolchain_search_sequence = [toolchain_gcc, toolchain_intel, toolchain_pgi]

	for candidate in toolchain_search_sequence:
		if conf.check_for_toolchain(candidate, "C++", "CXX", ".cpp"):
			toolchain = candidate
			break

	if not have_toolchain():
		print("error: could not identify toolchain on your system")
		Exit(1)

	if check_c and not conf.check_for_toolchain(toolchain, "C", "CC", ".c"):
		print("error: C toolchain does not match the identified C++ toolchain")
		Exit(1)

	myenv = conf.Finish()

	def add_flag_if_supported(env, tool, extension, flag):
		def check_flag_test(context, tool, extension, flag):
			## ...
		## ...

	def add_to_cflags_if_supported(env, flag):
		## ...
	
	def add_to_ccflags_if_supported(env, flag):
		## ...

	def add_to_cxxflags_if_supported(env, flag):
		## ...

	if using_gcc() or using_intel() or using_pgi():
		## ...
	
	if _has_option('c++11'):
		## ...

	conf = Configure(myenv)
	libdeps.setup_conftests(conf)

	return conf.Finish()


############### options and options

options = { }

# build output
add_option("mute" , "reduce screen noise", 0, False)
# installation/packaging
add_option("prefix", "installation prefix", 1, False, default=DEFAULT_INSTALL_DIR)
add_option("distname", "distribution name (1.0)", 1, False)
# linking options
add_option("release", "release build", 0, True)
add_option("static", "fully static build", 0, False)
# base compile flags
add_option("64", "whether to force 64 bit", 0, True, "force64")
add_option("32", "whether to force 32 bit", 0, True, "force32")
# compiler stuff
add_option("cxx", "compiler to use for c++", 1, True)
add_option("cc", "compiler to use for c", 1, True)
add_option("cc-use-shell-environment", "use $CC from shell for C compiler", 0, False)
add_option("cxx-use-shell-environment", "use $CXX from shell for C++ compiler", 0, False)
add_option("ld", "linker to use", 1, True)
add_option("c++11", "enable c++11 support (necessary)", 0, True)
add_option("cpppath", "Include path if you have headers in a nonstandard directory", 1, True)
add_option("libpath", "Library path if you have libraries in a nonstandard directory", 1, True)
add_option("extrapath", "comma separated list of additional paths, static linking", 1, True)
add_option("extrapathdyn", "comma separated list of additional paths, dynamic linking", 1, True)
add_option("extralib", "comma separated list of additional libraries", 1, True)
add_option("boost-compiler", "compiler used for boost (gcc46)", 1, True, "boostCompiler")
add_option("boost-version", "boost version for linking (1_49_0)", 1, True, "boostVersion")
# debug options
add_option("dbg", "debug build, no optimization", 0, True, "debugBuild")
add_option("dbg1", "debug build, no optimization, debug logging", 0, True, "debugBuildAndLog1")
add_option("use-papi", "Enable PAPI profiling", 0, False)
add_option("use-gpu", "Enable GPU support", 0, False)

printLocalInfo()

# do not build if help is called
if GetOption('help'):
	Return()

################ environment

boost_libs = ["system", "filesystem", "timer", "chrono"]
hdf5_libs = ["hdf5", "z"]
mpi_libs = ["mpi_cxx", "mpi"]
cuda_libs = ["cudart", "cufft", "cudadevrt", "nvToolsExt"]
tiff_libs = ["tiff"]
papi_libs = ["papi"]
other_libs = ["m"]

linux64 = True
force32 = _has_option("force32")
force64 = _has_option("force64")
if force32:
	print("32-bit is not supported. Assuming 64-bit")
if not force64 or not force32:
	print("Building for 64-bit")
msarch = "amd64"

release_build = _has_option("release")
using_debug = _has_option("debugBuild") or _has_option("debugBuildAndLog1")
if release_build and using_debug:
	print("error: release build cannot have debug enabled. make a choice first")
	Exit(1)

static = _has_option("static")

variant_dir = ".build"

# initialize the environment with gathered stuff so far
env = Environment(BUILD_DIR=variant_dir,
					DIST_ARCHIVE_SUFFIX = '.tar',
					EXTRAPATH = get_option("extrapath"),
					#PYTHON = utils.find_python(),
					TARGET_ARCH = msarch,
					PYSYSPLATFORM = os.sys.platform)

conf = Configure(env)

env['_LIBDEPS'] = '$_LIBDEPS_OBJS'

## the decider
env.Decider('MD5-timestamp')

if _has_option("mute"):
	env.Append(CCCOMSTR = "Compiling $TARGET")
	env.Append(CXXCOMSTR = env["CCCOMSTR"])
	env.Append(SHCCCOMSTR = env["CCCOMSTR"])
	env.Append(SHCXXCOMSTR = env["SHCCCOMSTR"])
	env.Append(LINKCOMSTR = "Linking $TARGET")
	env.Append(SHLINKCOMSTR = env["LINKCOMSTR"])
	env.Append(ARCOMSTR = "Generating library $TARGET")

boost_c = get_option("boostCompiler")
if boost_c is None:
	boost_c = ""
else:
	boost_c = '-' + boost_c

boost_v = get_option("boostVersion")
if boost_v is None:
	boost_v = ""
else:
	boost_v = '-' + boost_v

env['EXTRACPPPATH'] = [ ]
env['EXTRALIBPATH'] = [ ]

if _has_option("extrapath"):
	add_extra_paths(get_option("extrapath"))
if _has_option("extrapathdyn"):
	add_extra_paths(get_option("extrapathdyn"))
if _has_option("extralib"):
	add_extra_libs(get_option("extralib"))


class InstallSetup:
	binaries = False
	libraries = False
	headers = False

	def __init__(self):
		self.default()
	
	def default(self):
		self.binaries = True
		self.libraries = False
		self.headers = False

install_setup = InstallSetup()

if _has_option("full"):
	install_setup.headers = True
	install_setup.libraries = True


if "uname" in dir(os):
	processor = os.uname()[4]
else:
	processor = "x86_64"

env['PROCESSOR_ARCHITECTURE'] = processor

install_dir = DEFAULT_INSTALL_DIR
nix_lib_prefix = "lib"

if _has_option("prefix"):
	install_dir = get_option("prefix")

env.Append(LIBS = ['m'])
if os.uname()[4] == "x86_64":
	linux64 = True
	nix_lib_prefix = "lib64"
	env.Append(EXTRALIBPATH = ["/usr/lib64", "/lib64"])
	env.Append(LIBS = ["pthread"])

env.Append(CCFLAGS = ["-mavx -msse4.1 -msse4.2 -mssse3"])

if static:
	env.Append(LINKFLAGS = " -static ")

if using_debug:
	env.Append(CCFLAGS = ["-g -DDEBUG"])
else:
	env.Append(CCFLAGS = ["-O3 -DNDEBUG"])

if force64:
	env.Append(CCFLAGS = "-m64")
	env.Append(CCFLAGS = "-m64")

env.Append(CPPPATH = ['$EXTRACPPPATH'], LIBPATH = ['$EXTRALIBPATH'])

env = do_configure(env)
test_env = env.Clone()
test_env.Append(CPPPATH = ["../"])


if not env.GetOption('clean'):
	if not conf.CheckCXX():
		print('!! Your compiler and/or environment is not correctly configured.')
		Exit(0)

	if not conf.CheckFunc('printf'):
		print('!! Your compiler and/or environment is not correctly set.')
		Exit(0)

	if not conf.CheckCXXHeader('iostream'):
		print "You need 'iostream' to build"
		Exit(1)

	if not conf.CheckCXXHeader('vector'):
		print("You need 'vector' to build")
		Exit(1)

	if not conf.CheckLib('m'):
		print("You need 'libm' library to build")
		Exit(1)
	else:
		conf.env.Append('-lm')

	if not conf.CheckLib('libboost_system'):
		print("What you need is 'libboost_system'")
		conf.env['libboost_system'] = "no"
	else:
		conf.env.Append('-lboost_system')

	if 'CXX' in os.environ:
		conf.env.Replace(CXX = os.environ['CXX'])
		print(">> Using compiler " + os.environ['CXX'])

	if 'CXXFLAGS' in os.environ:
		conf.env.Append(CCFLAGS = os.environ['CXXFLAGS'])
		print(">> Appending custom build flags : " + os.environ['CXXFLAGS'])
						    
	if 'LDFLAGS' in os.environ:
		conf.env.Append(LINKFLAGS = os.environ['LDFLAGS'])
		print(">> Appending custom link flags : " + os.environ['LDFLAGS'])

	env = conf.Finish()


boost_dir = "/usr/local/boost_1_49_0"
hdf5_dir = "/usr/local/hdf5-1.8.9"
tiff_dir = "/usr/local/tiff-4.0.2"
mpi_dir = "/usr/local/openmpi"
papi_dir = "/usr/local/papi-5.1.0"

cuda_dir = "/usr/local/cuda"

