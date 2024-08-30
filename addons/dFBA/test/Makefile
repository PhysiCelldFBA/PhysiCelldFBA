CC := g++

# Check for environment definitions of compiler
ifdef PHYSICELL_CPP
	CC := $(PHYSICELL_CPP)
endif

ARCH := native # best auto-tuning

# For debugging with gdb just add: -g3 -ggdb -O3
CFLAGS      := -march=$(ARCH) -O3 -fomit-frame-pointer -mfpmath=both -fopenmp -m64 -std=c++11

# Libraries need: SBML and coin-Clp
DEPS_DIR    := ../ext
DEPS_CFLAGS := -I$(DEPS_DIR)/libsbml/include -I$(DEPS_DIR)/coin-or/include -I../src
LD_PATH     := -L$(DEPS_DIR)/libsbml/lib -L$(DEPS_DIR)/coin-or/lib/
# sudo apt install libbz2-dev
# sudo apt install libxml2-dev

## "linux-x64", "linux-x86", "win64", "win32", "osx"
PKG_ARCH := linux-x64

LIBS        := -lClp-static -lCoinUtils-static -llapack -lsbml-static -lxml2 -lbz2 -lz

ifeq ($(OS),Windows_NT)
	WIN_LD_FLAG := -LC:\msys64\mingw64\bin
	ARCH := win64
else
	WIN_LD_FLAG := 
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Darwin)
		UNAME_P := $(shell uname -p)
		var := $(shell which $(CC) | xargs file)
		ifeq ($(lastword $(var)),arm64)
		  CFLAGS := -march=$(ARCH) -O3 -fomit-frame-pointer -fopenmp -m64 -std=c++11
		endif
	endif
endif

COMPILE_COMMAND := $(CC) $(CFLAGS) $(DEPS_CFLAGS)

SRC_DIR := "./"

vpath %.o     $(SRC_DIR)
vpath %.cpp   $(SRC_DIR) 
vpath %.h     $(SRC_DIR) 

# ALL_OBJECTS := dfba_Metabolite.o dfba_Reaction.o dfba_Solution.o dfba_Model.o dfba_intracellular.o
ALL_OBJECTS := dfba_Metabolite.o dfba_Reaction.o dfba_Solution.o dfba_Model.o 

make_lib: libFBA $(ALL_OBJECTS)
all: make_lib test_sbml test_clp test_fba

# $(ALL_OBJECTS): %.o: %.cpp %.h
# 	$(COMPILE_COMMAND) $(LD_PATH) $(WIN_LD_FLAG) $(LIBS) -c $< -o $(BUILD_DIR)/$(@F)

test_sbml: libFBA $(ALL_OBJECTS) ./main_test_sbml.cpp
	$(COMPILE_COMMAND) -o test_sbml $(ALL_OBJECTS) main_test_sbml.cpp $(LD_PATH) $(WIN_LD_FLAG) $(LIBS)
	@echo created $(PROGRAM_NAME)
	@echo

test_clp: libFBA $(ALL_OBJECTS) ./main_test_clp.cpp
	$(COMPILE_COMMAND) -o test_clp $(ALL_OBJECTS) main_test_clp.cpp $(LD_PATH) $(WIN_LD_FLAG) $(LIBS)
	@echo created $(PROGRAM_NAME)
	@echo

test_fba: libFBA $(ALL_OBJECTS) ./main_test_fba.cpp
	$(COMPILE_COMMAND) -o test_fba $(ALL_OBJECTS) main_test_fba.cpp $(LD_PATH) $(WIN_LD_FLAG) $(LIBS)
	@echo created $(PROGRAM_NAME)
	@echo


# FBA addon modules
libFBA: 
	python3 ../../../beta/setup_fba.py --pkg libsbml --arch $(PKG_ARCH)
	python3 ../../../beta/setup_fba.py --pkg coin-or --arch $(PKG_ARCH)

dfba_intracellular.o: ../src/dfba_intracellular.cpp
	$(COMPILE_COMMAND) $(ADDFLAGS) $(INC) $(LIBS) -c ../src/dfba_intracellular.cpp -o ./dfba_intracellular.o

dfba_Model.o: ../src/dfba_Model.cpp
	$(COMPILE_COMMAND) $(ADDFLAGS) $(INC) $(LIBS) -c ../src/dfba_Model.cpp -o ./dfba_Model.o

dfba_Reaction.o: ../src/dfba_Reaction.cpp
	$(COMPILE_COMMAND) $(ADDFLAGS) $(INC) $(LIBS) -c ../src/dfba_Reaction.cpp -o ./dfba_Reaction.o

dfba_Metabolite.o: ../src/dfba_Metabolite.cpp
	$(COMPILE_COMMAND) $(ADDFLAGS) $(INC) $(LIBS) -c ../src/dfba_Metabolite.cpp -o ./dfba_Metabolite.o

dfba_Solution.o: ../src/dfba_Solution.cpp
	$(COMPILE_COMMAND) $(ADDFLAGS) $(INC) $(LIBS) -c ../src/dfba_Solution.cpp -o ./dfba_Solution.o


clean:
	rm -f *.o
	rm -f ./test_*
