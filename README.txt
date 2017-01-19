Dependencies:

	libhdf5 - change INCLUDE_DIRS and EXTERNAL_LIB_DIR in makefiles/Makefile.std to the location of the library and header files.
	gcc and OpenGL

Compilation:
	From the 'gaps' directory: make clean; make

This directory contains all code for the GAPS software library.
There are several subdirectories:

    pkgs - source and include files for all packages (software libraries).
    apps - source files for several application and example programs. 
    makefiles - unix-style make file definitions
    vc - visual studio solution files
    lib - archive library (.lib) files (created during compilation).
    bin - executable files (created during compilation).

The software is distributed under the MIT license (see LICENSE.txt)
and thus can be used for any purpose without warranty, any liability,
or any suport of any kind.



