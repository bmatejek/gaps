#
# Application name and list of source files.
#

NAME=skeleton
CCSRCS=$(NAME).cpp



#
# Dependency libraries
#

ifeq ($(shell uname), Linux)
PKG_LIBS=-lRNDataStructures -lR3Graphics -lR3Shapes -lR2Shapes -lRNBasics -lpng -ljpeg -lfglut -lX11 -lhdf5_cpp -lhdf5
else
PKG_LIBS=-lRNDataStructures -lR3Graphics -lR3Shapes -lR2Shapes -lRNBasics -lfglut -ljpeg -lpng -lhdf5_cpp -lhdf5
endif

#
# R3 application makefile
#

include ../../makefiles/Makefile.apps


