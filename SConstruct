env = Environment(CPPPATH = ['#/include'],
                  CC = 'g++',
                  CCFLAGS = '-O3 -std=c++17',
                  LIBS = [],
                  LIBPATH = [])

Export('env')

env.SConscript('src/SConscript', variant_dir="build", duplicate=0)
