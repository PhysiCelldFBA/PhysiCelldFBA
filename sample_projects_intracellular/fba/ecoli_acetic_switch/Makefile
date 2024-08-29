VERSION := $(shell grep . VERSION.txt | cut -f1 -d:)
PROGRAM_NAME := ecoli-dfba

CC := g++
# CC := g++-mp-7 # typical macports compiler name
# CC := g++-7 # typical homebrew compiler name

# Check for environment definitions of compiler
# e.g., on CC = g++-7 on OSX
ifdef PHYSICELL_CPP
	CC := $(PHYSICELL_CPP)
endif

ARCH := native # best auto-tuning
# ARCH := core2 # a reasonably safe default for most CPUs since 2007
# ARCH := corei7
# ARCH := corei7-avx # earlier i7
# ARCH := core-avx-i # i7 ivy bridge or newer
# ARCH := core-avx2 # i7 with Haswell or newer
# ARCH := nehalem
# ARCH := westmere
# ARCH := sandybridge # circa 2011
# ARCH := ivybridge   # circa 2012
# ARCH := haswell     # circa 2013
# ARCH := broadwell   # circa 2014
# ARCH := skylake     # circa 2015
# ARCH := bonnell
# ARCH := silvermont
# ARCH := skylake-avx512
# ARCH := nocona #64-bit pentium 4 or later


# For debugging with gdb just add: -g3 -ggdb -O3
# -DCLP_BUILD is an additional flags to link againt clp
# CFLAGS := -march=$(ARCH) -ggdb -fomit-frame-pointer -mfpmath=both -fopenmp -m64 -std=c++11 -DCLP_BUILD
CFLAGS := -march=$(ARCH) -fomit-frame-pointer -mfpmath=both -fopenmp -m64 -std=c++11 -DCLP_BUILD


# Libraries need: SBML and coin-Clp

CLP_FOLDER := ./addons/coin-clp

## "linux-x64", "linux-x86", "win64", "win32", "osx"
PKG_ARCH := linux-x64

# if we're on Windows, we need to provide a (agreed upon) path to some of the dynamic libs
ifeq ($(OS),Windows_NT)
	OSFLAG += -D WIN32
	OMP_LIB := 
#	LIBFBA_DIR := C:\Users\heiland\libroadrunner\roadrunner-win64-vs14-cp35m
#	LIBFBA_LIBS := C:\Users\heiland\libroadrunner\roadrunner-win64-vs14-cp35m/bin
#	LIBFBA_CFLAGS := -I${LIBFBA_DIR}/include/rr/C
	CFLAGS := -march=$(ARCH) -fomit-frame-pointer -fopenmp -m64 -std=c++11 -D ADDON_PHYSIDFBA
	ifeq ($(PROCESSOR_ARCHITECTURE),AMD64)
		OSFLAG += -D AMD64
	endif
	ifeq ($(PROCESSOR_ARCHITECTURE),x86)
		OSFLAG += -D IA32
	endif
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

COMPILE_COMMAND := $(CC) $(CFLAGS) $(LIBFBA_CFLAGS)

BioFVM_OBJECTS := BioFVM_vector.o BioFVM_mesh.o BioFVM_microenvironment.o BioFVM_solvers.o BioFVM_matlab.o \
BioFVM_utilities.o BioFVM_basic_agent.o BioFVM_MultiCellDS.o BioFVM_agent_container.o

PhysiCell_core_OBJECTS := PhysiCell_phenotype.o PhysiCell_cell_container.o PhysiCell_standard_models.o \
PhysiCell_cell.o PhysiCell_custom.o PhysiCell_utilities.o PhysiCell_constants.o PhysiCell_basic_signaling.o \
PhysiCell_signal_behavior.o PhysiCell_rules.o

PhysiCell_module_OBJECTS := PhysiCell_SVG.o PhysiCell_pathology.o PhysiCell_MultiCellDS.o PhysiCell_various_outputs.o \
PhysiCell_pugixml.o PhysiCell_settings.o  PhysiCell_geometry.o

PhysiCelldFBA_addon_OBJECTS := dfba_Metabolite.o dfba_Reaction.o dfba_Solution.o dfba_Model.o dfba_intracellular.o


# put your custom objects here (they should be in the custom_modules directory)

PhysiCell_custom_module_OBJECTS := cutom.o

pugixml_OBJECTS := pugixml.o

PhysiCell_OBJECTS := $(BioFVM_OBJECTS)  $(pugixml_OBJECTS) $(PhysiCell_core_OBJECTS) $(PhysiCell_module_OBJECTS)
ALL_OBJECTS := $(PhysiCell_OBJECTS) $(PhysiCelldFBA_addon_OBJECTS) $(PhysiCell_custom_module_OBJECTS)


# compile the project

# all: libFBA main.cpp $(ALL_OBJECTS)

all: main.cpp $(ALL_OBJECTS)
	$(COMPILE_COMMAND) $(INC) -o $(PROGRAM_NAME) $(ALL_OBJECTS) main.cpp $(LIBS)
	make name 

name:
	@echo ""
	@echo "Executable name is" $(PROGRAM_NAME)
	@echo ""

# PhysiCell core components

PhysiCell_phenotype.o: ./core/PhysiCell_phenotype.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_phenotype.cpp

PhysiCell_digital_cell_line.o: ./core/PhysiCell_digital_cell_line.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_digital_cell_line.cpp

PhysiCell_cell.o: ./core/PhysiCell_cell.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_cell.cpp

PhysiCell_cell_container.o: ./core/PhysiCell_cell_container.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_cell_container.cpp

PhysiCell_standard_models.o: ./core/PhysiCell_standard_models.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_standard_models.cpp

PhysiCell_utilities.o: ./core/PhysiCell_utilities.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_utilities.cpp

PhysiCell_custom.o: ./core/PhysiCell_custom.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_custom.cpp
	
PhysiCell_constants.o: ./core/PhysiCell_constants.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_constants.cpp 

PhysiCell_signal_behavior.o: ./core/PhysiCell_signal_behavior.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_signal_behavior.cpp 

PhysiCell_basic_signaling.o: ./core/PhysiCell_basic_signaling.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_basic_signaling.cpp
	
PhysiCell_geometry.o: ./modules/PhysiCell_geometry.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_geometry.cpp 

PhysiCell_rules.o: ./core/PhysiCell_rules.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_rules.cpp 

# BioFVM core components (needed by PhysiCell)

BioFVM_vector.o: ./BioFVM/BioFVM_vector.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_vector.cpp

BioFVM_agent_container.o: ./BioFVM/BioFVM_agent_container.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_agent_container.cpp

BioFVM_mesh.o: ./BioFVM/BioFVM_mesh.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_mesh.cpp

BioFVM_microenvironment.o: ./BioFVM/BioFVM_microenvironment.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_microenvironment.cpp

BioFVM_solvers.o: ./BioFVM/BioFVM_solvers.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_solvers.cpp

BioFVM_utilities.o: ./BioFVM/BioFVM_utilities.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_utilities.cpp

BioFVM_basic_agent.o: ./BioFVM/BioFVM_basic_agent.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_basic_agent.cpp

BioFVM_matlab.o: ./BioFVM/BioFVM_matlab.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_matlab.cpp

BioFVM_MultiCellDS.o: ./BioFVM/BioFVM_MultiCellDS.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_MultiCellDS.cpp

pugixml.o: ./BioFVM/pugixml.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/pugixml.cpp

# standard PhysiCell modules

PhysiCell_SVG.o: ./modules/PhysiCell_SVG.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_SVG.cpp

PhysiCell_pathology.o: ./modules/PhysiCell_pathology.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_pathology.cpp

PhysiCell_MultiCellDS.o: ./modules/PhysiCell_MultiCellDS.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_MultiCellDS.cpp

PhysiCell_various_outputs.o: ./modules/PhysiCell_various_outputs.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_various_outputs.cpp

PhysiCell_pugixml.o: ./modules/PhysiCell_pugixml.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_pugixml.cpp

PhysiCell_settings.o: ./modules/PhysiCell_settings.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_settings.cpp

# user-defined PhysiCell modules

# dFBA addon modules

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

# user-defined PhysiCell modules
cutom.o: ./custom_modules/cutom.cpp
	$(COMPILE_COMMAND) -c ./custom_modules/cutom.cpp


# cleanup

reset:
	rm -f *.cpp
	cp ./sample_projects/Makefile-default Makefile
	rm -f ./custom_modules/*
	touch ./custom_modules/empty.txt
	touch ALL_CITATIONS.txt
	touch ./core/PhysiCell_cell.cpp
	rm ALL_CITATIONS.txt
	cp ./config/PhysiCell_settings-backup.xml ./config/PhysiCell_settings.xml
	touch ./config/empty.txt

clean:
	rm -f *.o
	rm -f $(PROGRAM_NAME)*

data-cleanup:
	rm -f *.mat
	rm -f *.xml
	rm -f *.svg
	rm -f ./output/*
	touch ./config/empty.txt

# archival

checkpoint: 
	zip -r $$(date +%b_%d_%Y_%H%M).zip Makefile *.cpp *.h config/*.xml custom_modules/* 

zip:
	zip -r latest.zip Makefile* *.cpp *.h BioFVM/* config/* core/* custom_modules/* matlab/* modules/* sample_projects/*
	cp latest.zip $$(date +%b_%d_%Y_%H%M).zip
	cp latest.zip VERSION_$(VERSION).zip
	mv *.zip archives/

tar:
	tar --ignore-failed-read -czf latest.tar Makefile* *.cpp *.h BioFVM/* config/* core/* custom_modules/* matlab/* modules/* sample_projects/*
	cp latest.tar $$(date +%b_%d_%Y_%H%M).tar
	cp latest.tar VERSION_$(VERSION).tar
	mv *.tar archives/

unzip:
	cp ./archives/latest.zip .
	unzip latest.zip

untar:
	cp ./archives/latest.tar .
	tar -xzf latest.tar
