#
# Application name 
#

NAME=skel_dist



#
# Source files
#

CCSRCS=$(NAME).cpp

#
# Libraries
#

ifeq ($(shell uname), Linux)
PKG_LIBS=-lRNDataStructures -lR3Graphics -lR3Shapes -lR2Shapes -lRNBasics -lpng -ljpeg -lhdf5_cpp -lhdf5
else
PKG_LIBS=-lRNDataStructures -lR3Graphics -lR3Shapes -lR2Shapes -lRNBasics -lpng -ljpeg -lhdf5_cpp -lhdf5
endif


#
# Include standard makefile
#

include ../../makefiles/Makefile.apps















