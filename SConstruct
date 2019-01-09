import os, sys
from os.path import join



##### define functions
def uniqueCheckLib(conf, lib):
   """ Checks for a library and appends it to env if not already appended. """
   if conf.CheckLib(lib, autoadd=0, language="C++"):
      conf.env.AppendUnique(LIBS = [lib])
      return True
   else:
      print("ERROR: Library '" + lib + "' not found!")
      Exit(1)
#######



# add possibility to add a debug visu
vars = Variables('custom.py')
vars.Add(BoolVariable('visu', 'Set to 1 for enabling debug visu', 0))

env = Environment(variables = vars, ENV = os.environ)
conf = Configure(env)


# ====== precice ======

preciceRoot = os.getenv ('PRECICE_ROOT')

if preciceRoot:
    print("PRECICE_ROOT defined, preCICE was probably build from source")
    print('Using environment variable PRECICE_ROOT = ' + preciceRoot)
    env.Append(CPPPATH = [os.path.join(preciceRoot, 'src')])
    env.Append(LIBPATH = [os.path.join(preciceRoot, 'build/last')]) 
    env.Append(CPPDEFINES = ['PRECICE_USE_MPI'])
    uniqueCheckLib(conf, "precice") 
else:
    print("PRECICE_ROOT not defined. Using pkg_config to find libprecice.")
    try:
        uniqueCheckLib(conf, "precice")
    except Exception():
        print("Did you forget to define PRECICE_ROOT?")
        Exit(-1)


# do debug build?
debug = ARGUMENTS.get('debug', 0)

# parallel
env.Replace(CXX='mpic++')

# serial
# env.Replace(CXX='g++')

# define some general compiler flags
env.Append(
    CXXFLAGS=[
        "-Wall",
        "-Wextra",
        "-pedantic",
        "-std=c++11",
    ],
    LIBS=[
        "mpi",
        "SDL2",
    ]
)

# add flags for debug and release build
if debug == 0:
    env['CXXFLAGS'] += ["-O3"]
else:
    env['CXXFLAGS'] += [
        "-g3",
        "-O0",
    ]

# call SConscript to actually build the project after setting up the environment
env.SConscript("./FluidSolver/SConscript", exports='env', variant_dir='./build/FluidSolver', duplicate=0)
env.SConscript("./Magrathea/SConscript", exports='env', variant_dir='./build/Magrathea', duplicate=0)
