ifndef APPLYFORCE
APPLYFORCE=APPLYFORCE

ifndef CLEANFOLDER
CLEANFOLDER=APPLYFORCE
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the object files to be compiled for this utility
APPLYFORCE_OBJECTS=ApplyForce.o   #initGraphics.o

# the libraries this utility depends on
APPLYFORCE_LIBS=sceneObject integratorSparse integrator elasticForceModel forceModel sparseMatrix loadList insertRows lighting configFile volumetricMesh getopts camera graph isotropicHyperelasticFEM stvk corotationalLinearFEM polarDecomposition minivector matrixIO massSpringSystem objMesh imageIO sparseSolver uniqueIntegerID 

# the headers in this library
#APPLYFORCE_HEADERS=initGraphics.h

APPLYFORCE_LINK=$(addprefix -l, $(APPLYFORCE_LIBS)) $(SPOOLES_LIB) $(PARDISO_LIB) $(GLUI_LIB) $(IMAGE_LIBS) $(STANDARD_LIBS)

APPLYFORCE_OBJECTS_FILENAMES=$(addprefix $(R)/utilities/ApplyForce/, $(APPLYFORCE_OBJECTS))
APPLYFORCE_HEADER_FILENAMES=$(addprefix $(R)/utilities/ApplyForce/, $(APPLYFORCE_HEADERS))
APPLYFORCE_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(APPLYFORCE_LIBS))
APPLYFORCE_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(APPLYFORCE_LIBS))

include $(APPLYFORCE_LIB_MAKEFILES) ../../libraries/glui/Makefile

all: $(R)/utilities/ApplyForce/ApplyForce

CURRENT_DIR = $(shell pwd)
$(R)/utilities/ApplyForce/ApplyForce: $(APPLYFORCE_OBJECTS_FILENAMES)
	$(CXXLD) $(LDFLAGS) $(APPLYFORCE_OBJECTS) $(APPLYFORCE_LINK) -Wl,-rpath,$(CURRENT_DIR)/$(GLUI_DIR)/lib -o $@; cp $@ $(R)/utilities/bin/

$(APPLYFORCE_OBJECTS_FILENAMES): %.o: %.cpp $(APPLYFORCE_LIB_FILENAMES) $(APPLYFORCE_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $(GLUI_INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), APPLYFORCE)
clean: cleanApplyForce
endif

deepclean: cleanApplyForce

cleanApplyForce:
	$(RM) $(APPLYFORCE_OBJECTS_FILENAMES) $(R)/utilities/ApplyForce/ApplyForce

endif

