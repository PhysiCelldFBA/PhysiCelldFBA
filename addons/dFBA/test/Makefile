# Compiler
CC := g++

ifdef PHYSICELL_CPP
	CC := $(PHYSICELL_CPP)
endif

# Architecture
ARCH := native

# Flags
CFLAGS := -march=$(ARCH) -O3 -fomit-frame-pointer -fopenmp -m64 -std=c++11

# Paths
CUR_DIR := $(shell pwd)/..
DEPS_DIR := $(CUR_DIR)/ext
SRC_DIR := ../src

INCLUDES := -I$(DEPS_DIR)/libsbml/include -I$(DEPS_DIR)/coin-or/include -I$(SRC_DIR)

LIBDIRS := -L$(DEPS_DIR)/libsbml/lib -L$(DEPS_DIR)/coin-or/lib

# LIBS := $(DEPS_DIR)/coin-or/lib/libClp.a $(DEPS_DIR)/coin-or/lib/libCoinUtils.a -llapack -lsbml-static -lxml2 -lbz2 -lz
LIBS := \
    $(DEPS_DIR)/coin-or/lib/libOsiClp.a \
    $(DEPS_DIR)/coin-or/lib/libOsi.a \
    $(DEPS_DIR)/coin-or/lib/libClp.a \
    $(DEPS_DIR)/coin-or/lib/libCoinUtils.a \
	$(DEPS_DIR)/libsbml/lib/libsbml-static.a \
    -llapack -lxml2 -lbz2 -lz


# Platform-specific tweaks
ifeq ($(OS),Windows_NT)
	WIN_LD_FLAG := -LC:\msys64\mingw64\bin
else
	WIN_LD_FLAG :=
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Darwin)
		var := $(shell which $(CC) | xargs file)
		ifeq ($(lastword $(var)),arm64)
			CFLAGS := -march=$(ARCH) -O3 -fomit-frame-pointer -fopenmp -std=c++11
		endif
	endif
endif

# Build
COMPILE_COMMAND := $(CC) $(CFLAGS) $(INCLUDES)

# Objects
OBJECTS := dfba_Metabolite.o dfba_Reaction.o dfba_Solution.o dfba_Model.o

# Targets
all: libFBA test_sbml test_clp test_fba

libFBA:
	python3 ../../../beta/setup_fba.py --pkg libsbml --arch linux-x64
	python3 ../../../beta/setup_fba.py --pkg coin-or --arch linux-x64

%.o: $(SRC_DIR)/%.cpp
	$(COMPILE_COMMAND) -c $< -o $@

test_sbml: libFBA $(OBJECTS) main_test_sbml.cpp
	$(COMPILE_COMMAND) -o $@ $(OBJECTS) main_test_sbml.cpp $(LIBDIRS) $(WIN_LD_FLAG) $(LIBS)
	@echo "Created $@"

test_clp: libFBA $(OBJECTS) main_test_clp.cpp
	$(COMPILE_COMMAND) -o $@ $(OBJECTS) main_test_clp.cpp $(LIBDIRS) $(WIN_LD_FLAG) $(LIBS)
	@echo "Created $@"

test_warmstart: libFBA $(ALL_OBJECTS) ./main_test_warmstart.cpp
	$(COMPILE_COMMAND) -o test_warmstart $(ALL_OBJECTS) main_test_warmstart.cpp $(LD_PATH) $(WIN_LD_FLAG) $(LIBS)
	@echo "created test_warmstart"
	@echo

test_fba: libFBA $(OBJECTS) main_test_fba.cpp
	$(COMPILE_COMMAND) -o $@ $(OBJECTS) main_test_fba.cpp $(LIBDIRS) $(WIN_LD_FLAG) $(LIBS)
	@echo "Created $@"

clean:
	rm -f *.o test_sbml test_clp test_fba
