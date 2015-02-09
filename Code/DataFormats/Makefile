include ../make.inc
#################################################
### Defining Compiling variables              ###
#################################################



CXX           = g++ -Wno-deprecated -Wall
LD            = g++ -Wno-deprecated -Wall
INSTALL	      = /usr/bin/install
SHELL = /bin/sh

####################################
###                              ###
####################################

override CXXFLAGS += -I$(ROOTSYS)/include -I./ -c
LDFLAGS       =  
ROOTLIBS      =  -L$(ROOTSYS)/lib -L/usr/lib/ -L/lib/i686/  -lCore -lCint -lHist -lGraf  -lGraf3d -lGpad -lTree -lRint -lReflexDict -lReflex -lPostscript -lMatrix -lPhysics -lMinuit2 -lGui -LObj -lThread -rdynamic -Wl,--rpath $(ROOTSYS)/lib

DEFS  	      = -DSTANDALONE=1
LIBS          = $(ROOTLIBS) $(SHAREDLIBFLAGS)

SRCS    = $(wildcard *.cxx)
HDRS    = $(SRCS:.cxx=.h)

OBJS          = $(SRCS:.cxx=.o) 

LIB     = lib/libDataFormats.so


##########################################################
###                                                    ###
##########################################################

ROOTCINT: $(HDRS)  Linkdef.h
	 rootcint -f MyDict.cxx -c  $^

$(LIB): MyDict.cxx 
	g++ -shared -o $@ `root-config --ldflags` $(CXXFLAGS) $(LDFLAGS) $(ROOTLIBS) $^

.PHONY: clean install 

install: ROOTCINT $(LIB)

clean:
	@rm MyDict.*
	@rm lib/*.so