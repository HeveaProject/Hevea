vpath %.cpp ../src/

prefix = $(CURDIR)/..
blitzprefix = /usr/local/blitz-0.9
# Path where all the packages for Hevea are installed
HEVEA = $(prefix)
# Path where Blitz++ is installed
BLITZ = $(blitzprefix)

LIBSDIR = $(BLITZ)/lib

BZDIR = $(BLITZ)/include
MATDIR = $(HEVEA)/matrices/include
ARGDIR = $(HEVEA)/argstream/v0r1m0/include
UTIDIR = $(HEVEA)/utils/include
TORDIR = $(HEVEA)/torus/v0r1m2/include
ISODIR = $(HEVEA)/iso/v0r1m1/include
INTDIR = $(HEVEA)/integration/include
CYLDIR = $(HEVEA)/cylinder/v0r0m1/include
LOCDIR = $(HEVEA)/local/v0r1m0/include

INCLUDES = -I$(BZDIR) -I$(MATDIR) -I$(ARGDIR) -I$(TORDIR) -I$(UTIDIR) -I$(ISODIR) -I$(INTDIR) -I$(CYLDIR) -I$(LOCDIR)

CXX = g++

# Set DEBUGMODE to yes if debug compilation is desired
DEBUGMODE = no
# Set OPENMPMODE to no to disable the openMP parallel programming library 
OPENMPMODE = yes

ifeq ($(OPENMPMODE),yes)
OPENMPFLAG = -fopenmp
else
OPENMPFLAG = 
endif

ifeq ($(DEBUGMODE),yes)
# Flags for debugging
CXXFLAGS = $(OPENMPFLAG) -ftemplate-depth-30 -g -DBZ_DEBUG 
LDFLAGS = -g $(OPENMPFLAG)
else
# Flags for optimized executables. 
CXXFLAGS = -O3 $(OPENMPFLAG) -ftemplate-depth-30 
LDFLAGS = $(OPENMPFLAG)
endif

LIBS = -L$(LIBSDIR) -lblitz -lm

TARGETS = flattorus

all: $(TARGETS) 

%.o: %.cpp 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $^

$(TARGETS): $(TARGETS).o
	$(CXX) $(LDFLAGS) $@.o -o $@ $(LIBS)

clean:
	-rm -f *.o $(TARGETS)

