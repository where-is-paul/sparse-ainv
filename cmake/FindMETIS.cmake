# This module defines the following variables:
#
# METIS_INCLUDE_DIR 
# METIS_LIBRARIES   - list of libraries to link against when using METIS.
# METIS_FOUND       - True if METIS was found.

SET ( METIS_HOME "./extern/metis/" )

FIND_PATH( METIS_INCLUDE_DIR metis.h
    ${METIS_HOME}/include/
    DOC "The directory where METIS resides")

FIND_LIBRARY( METIS_LIBRARY
    NAME metis
    PATHS
    ${METIS_HOME}/lib    
    DOC "The METIS library")

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(METIS DEFAULT_MSG
    METIS_INCLUDE_DIR 
    METIS_LIBRARY )

IF( METIS_FOUND )
  set( METIS_INCLUDE_DIRS ${METIS_INCLUDE_DIR} )
  set( METIS_LIBRARIES ${METIS_LIBRARY} )
endif( METIS_FOUND )