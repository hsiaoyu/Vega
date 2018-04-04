ifndef DISPLACEMENT
DISPLACEMENT=DISPLACEMENT

ifndef CLEANFOLDER
CLEANFOLDER=DISPLACEMENT
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the object files to be compiled for this utility
DISPLACEMENT_OBJECTS=DisplaceBoundary.o   #initGraphics.o

# the libraries this utility depends on
DISPLACEMENT_LIBS=sceneObject integratorSparse integrator elasticForceModel forceModel sparseMatrix loadList insertRows lighting configFile volumetricMesh getopts camera graph isotropicHyperelasticFEM stvk corotationalLinearFEM polarDecomposition minivector matrixIO massSpringSystem objMesh imageIO sparseSolver uniqueIntegerID 

# the headers in this library
#DISPLACEMENT_HEADERS=initGraphics.h

DISPLACEMENT_LINK=$(addprefix -l, $(DISPLACEMENT_LIBS)) $(SPOOLES_LIB) $(PARDISO_LIB) $(GLUI_LIB) $(IMAGE_LIBS) $(STANDARD_LIBS)

DISPLACEMENT_OBJECTS_FILENAMES=$(addprefix $(R)/utilities/DisplaceBoundary/, $(DISPLACEMENT_OBJECTS))
DISPLACEMENT_HEADER_FILENAMES=$(addprefix $(R)/utilities/DisplaceBoundary/, $(DISPLACEMENT_HEADERS))
DISPLACEMENT_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(DISPLACEMENT_LIBS))
DISPLACEMENT_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(DISPLACEMENT_LIBS))

include $(DISPLACEMENT_LIB_MAKEFILES) ../../libraries/glui/Makefile

all: $(R)/utilities/DisplaceBoundary/DisplaceBoundary

CURRENT_DIR = $(shell pwd)
$(R)/utilities/DisplaceBoundary/DisplaceBoundary: $(DISPLACEMENT_OBJECTS_FILENAMES)
	$(CXXLD) $(LDFLAGS) $(DISPLACEMENT_OBJECTS) $(DISPLACEMENT_LINK) -Wl,-rpath,$(CURRENT_DIR)/$(GLUI_DIR)/lib -o $@; cp $@ $(R)/utilities/bin/

$(DISPLACEMENT_OBJECTS_FILENAMES): %.o: %.cpp $(DISPLACEMENT_LIB_FILENAMES) $(DISPLACEMENT_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $(GLUI_INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), DISPLACEMENT)
clean: cleanDisplaceBoundary
endif

deepclean: cleanDisplaceBoundary

cleanDisplaceBoundary:
	$(RM) $(DISPLACEMENT_OBJECTS_FILENAMES) $(R)/utilities/DisplaceBoundary/DisplaceBoundary

endif

