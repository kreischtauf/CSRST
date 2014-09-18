#
# Find Boost includes and library
#
# Boost
# It can be found at:
#     http://boost.org
#
# BOOST_INCLUDE_DIR - where to find hdf5.h
# BOOST_LIBRARY     - qualified libraries to link against.
# BOOST_FOUND       - do not attempt to use if "no" or undefined.

#SET(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a .so ${CMAKE_FIND_LIBRARY_SUFFIXES})
FIND_PATH(BOOST_INCLUDE_DIR boost/regex.hpp
  PATHS $ENV{BOOST_INCLUDE_PATH}
  NO_DEFAULT_PATH
)

FIND_PATH(BOOST_INCLUDE_DIR boost/regex.hpp
  /usr/local/include
  /usr/include
)

FIND_LIBRARY(BOOST_REGEX_LIBRARY boost_regex
  PATHS $ENV{BOOST_LIBRARY_PATH}
  NO_DEFAULT_PATH
)

FIND_LIBRARY(BOOST_REGEX_LIBRARY boost_regex
  /usr/lib
  /usr/local/lib
)

FIND_LIBRARY(BOOST_SERIALIZATION_LIBRARY boost_serialization
  PATHS $ENV{BOOST_LIBRARY_PATH}
  NO_DEFAULT_PATH
)

FIND_LIBRARY(BOOST_SERIALIZATION_LIBRARY boost_serialization
  /usr/lib
  /usr/local/lib
)

FIND_LIBRARY(BOOST_MPI_LIBRARY boost_mpi
  PATHS $ENV{BOOST_LIBRARY_PATH}
  NO_DEFAULT_PATH
)

FIND_LIBRARY(BOOST_MPI_LIBRARY boost_mpi
  /usr/lib
  /usr/local/lib
)

FIND_LIBRARY(BOOST_TIMER_LIBRARY boost_timer
  PATHS $ENV{BOOST_LIBRARY_PATH}
  NO_DEFAULT_PATH
)

FIND_LIBRARY(BOOST_TIMER_LIBRARY boost_timer
  /usr/lib
  /usr/local/lib
)

FIND_LIBRARY(BOOST_SYSTEM_LIBRARY boost_system
  PATHS $ENV{BOOST_LIBRARY_PATH}
  NO_DEFAULT_PATH
)

FIND_LIBRARY(BOOST_SYSTEM_LIBRARY boost_system
  /usr/lib
  /usr/local/lib
)

IF(BOOST_INCLUDE_DIR AND BOOST_REGEX_LIBRARY AND BOOST_SERIALIZATION_LIBRARY AND BOOST_MPI_LIBRARY AND BOOST_SYSTEM_LIBRARY AND BOOST_TIMER_LIBRARY)
    SET( BOOST_FOUND "YES" )
ENDIF(BOOST_INCLUDE_DIR AND BOOST_REGEX_LIBRARY AND BOOST_SERIALIZATION_LIBRARY AND BOOST_MPI_LIBRARY AND BOOST_SYSTEM_LIBRARY AND BOOST_TIMER_LIBRARY)

IF (BOOST_FOUND)
   IF (NOT BOOST_FIND_QUIETLY)
      MESSAGE(STATUS "Found BOOST: ${BOOST_REGEX_LIBRARY}; ${BOOST_SERIALIZATION_LIBRARY}; ${BOOST_MPI_LIBRARY}; ${BOOST_TIMER_LIBRARY}; ${BOOST_SYSTEM_LIBRARY}")
   ENDIF (NOT BOOST_FIND_QUIETLY)
ELSE (BOOST_FOUND)
   IF (BOOST_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find BOOST!")
   ENDIF (BOOST_FIND_REQUIRED)
ENDIF (BOOST_FOUND)
