
CXXFLAGS += -I${CLHEP_INCLUDE_DIR}

LDFLAGS += -L${CLHEP_LIB_DIR}

ROOTCFLAGS    = $(shell root-config --cflags --glibs )
ROOTLIBS      = $(shell root-config --libs)

%:%.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(ROOTCFLAGS) $(ROOTLIBS) -lRooFit -lHtml -lMinuit $< -o $@  

%TARGETS = kalFit
TARGETS = kalFit

all: $(TARGETS) 

clean:
	rm -rf $(TARGETS)

 NugenDeuteriumGen :  NugenDeuteriumGen.o
	g++ -o  NugenDeuteriumGen  NugenDeuteriumGen.o -L${ROOTSYS}/lib $(ROOTLIBS) $(ROOTCFLAGS) -lm -lc
 NugenDeuteriumGen.o :  NugenDeuteriumGen.cc
	g++ -c ${ROOTCFLAGS} -I${DOGS_PATH}/DCNuGen2 -I${DOGS_PATH}/DCBase -lRooFit -lHtml -lMinuit -I${DOGS_PATH}/DCGeo -I${DOGS_PATH}/DCValidity -I${DOGS_PATH}/DCEvent  NugenDeuteriumGen.cc

