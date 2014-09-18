#
# Find LAPACK includes and library
#
# LAPACK 
# It can be found at:
#     http://www.netlib.org/lapack/
#
# LAPACK_LIBRARY     - qualified libraries to link against.
# LAPACK_FOUND       - do not attempt to use if "no" or undefined.

SET(CMAKE_FIND_LIBRARY_SUFFIXES .lib .so .a ${CMAKE_FIND_LIBRARY_SUFFIXES})

FIND_LIBRARY(LAPACK_LIBRARY lapack
   /usr/local/lib
   /usr/lib
   $ENV{LAPACK_LIBRARY_PATH}
)

IF(LAPACK_LIBRARY)
    SET( LAPACK_FOUND "YES" )
ENDIF(LAPACK_LIBRARY)

IF (LAPACK_FOUND)
   IF (NOT LAPACK_FIND_QUIETLY)
      MESSAGE(STATUS "Found LAPACK: ${LAPACK_LIBRARY}")
   ENDIF (NOT LAPACK_FIND_QUIETLY)
ELSE (LAPACK_FOUND)
   IF (LAPACK_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find LAPACK!")
   ENDIF (LAPACK_FIND_REQUIRED)
ENDIF (LAPACK_FOUND)
