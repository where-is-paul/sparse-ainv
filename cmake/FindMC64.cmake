# This module defines the following variables:
#
# MC64_INCLUDE_DIR 
# MC64_LIBRARIES   - list of libraries to link against when using MC64.
# MC64_FOUND       - True if MC64 was found.

SET ( MC64_HOME "./extern/hsl_mc64/" )

FIND_PATH( MC64_INCLUDE_DIR hsl_mc64d.h
    ${MC64_HOME}/include/
    DOC "The directory where MC64 resides")

FIND_LIBRARY( MC64_LIBRARY
    NAME hsl_mc64
    PATHS
    ${MC64_HOME}/lib    
    DOC "The MC64 library")

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MC64 DEFAULT_MSG
    MC64_INCLUDE_DIR 
    MC64_LIBRARY )

IF( MC64_FOUND )
  set( MC64_INCLUDE_DIRS ${MC64_INCLUDE_DIR} )
  set( MC64_LIBRARIES ${MC64_LIBRARY} )
endif( MC64_FOUND )