#
# Find IPPL includes and library

# It can be found at:
#     http://amas.web.psi.ch/tools/IPPL/index.html
#
# Metis_LIBRARY     - qualified libraries to link against.
# Metis_FOUND       - do not attempt to use if "no" or undefined.

FIND_LIBRARY(Metis_LIBRARY metis
  /usr/local/lib
  /usr/lib
  $ENV{Metis_ROOT}/lib
)

IF(Metis_LIBRARY)
    SET( Metis_FOUND "YES" )
ENDIF(Metis_LIBRARY)

IF (Metis_FOUND)
   IF (NOT Metis_FIND_QUIETLY)
      MESSAGE(STATUS "Found Metis: ${Metis_LIBRARY}")
   ENDIF (NOT Metis_FIND_QUIETLY)
ELSE (Metis_FOUND)
   IF (Metis_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find Metis!")
   ENDIF (Metis_FIND_REQUIRED)
ENDIF (Metis_FOUND)
