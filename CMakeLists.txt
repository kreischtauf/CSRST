cmake_minimum_required (VERSION 2.6)

project(CSRST)

SET (CMAKE_CXX_FLAGS_RELEASE "-O3")
SET (CMAKE_CXX_FLAGS_DEBUG
     "-g -O2" CACHE STRING "Debug options." FORCE)
SET (CMAKE_CXX_COMPILER ${MPI_COMPILER})

SET (${PROJECT_NAME}_VERSION_MAJOR 0)
SET (${PROJECT_NAME}_VERSION_MINOR 1)
IF(EXISTS ".svn")
    include(FindSubversion)
    Subversion_WC_INFO(${CMAKE_CURRENT_SOURCE_DIR} OPAL)
ELSE(EXISTS ".svn")
    SET (OPAL_WC_REVISION 0)
    EXECUTE_PROCESS(
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMAND git rev-parse HEAD
        OUTPUT_VARIABLE GIT_OUT OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    SET (${PROJECT_NAME}_VERSION_PATCH \"${GIT_OUT}\")
ENDIF(EXISTS ".svn")

SET (PACKAGE_NAME \"CSRST\")
STRING (TOLOWER ${PACKAGE_NAME} PACKAGE)
SET (PACKAGE_TARNAME ${PACKAGE})
SET (PACKAGE_VERSION \"${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}\")
SET (VERSION ${PACKAGE_VERSION})
SET (PACKAGE_STRING ${PACKAGE_NAME}\ ${PACKAGE_VERSION})

SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

MACRO (MACRO_OPTIONAL_FIND_PACKAGE _name )
   OPTION(WITH_${_name} "Search for ${_name} package" OFF)
   if (WITH_${_name})
      FIND_PACKAGE(${_name} ${ARGN})
      set(DEFINE_WITH_${_name} "#define ${_name}_FOUND")
   else (WITH_${_name})
      set(${_name}_FOUND)
      set(${_name}_INCLUDE_DIR)
      set(${_name}_INCLUDES)
      set(${_name}_LIBRARY)
      set(${_name}_LIBRARIES)
   endif (WITH_${_name})
ENDMACRO (MACRO_OPTIONAL_FIND_PACKAGE)

find_package (IPPL REQUIRED
               PATHS /usr/local/include /usr/include $ENV{IPPL_PREFIX})
IF(${IPPL_LIBRARY} STREQUAL "ippl")
  MESSAGE(STATUS "Found IPPL: ${IPPL_LIBRARY_DIR}/libippl.a")
  MESSAGE(STATUS "            ${IPPL_INCLUDE_DIR}")
ENDIF(${IPPL_LIBRARY} STREQUAL "ippl")

find_package (H5Hut REQUIRED)
find_package (HDF5 REQUIRED)
find_package (TCLAP REQUIRED)
find_package (Boost REQUIRED)
find_package (Metis REQUIRED)
find_package (Lapack REQUIRED)
find_package (GSL REQUIRED)
macro_optional_find_package (OpenSSL)
#find_package (OpenSSL REQUIRED)
#find_package (SZIP REQUIRED)
find_package (Trilinos PATHS $ENV{TRILINOS_INCLUDE_PATH})

IF(NOT Trilinos_FOUND)
    MESSAGE (FATAL_ERROR "Could not find Trilinos!")
ENDIF()
# Echo trilinos build info just for fun
MESSAGE("\nFound Trilinos!  Here are the details: ")
MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
MESSAGE("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
MESSAGE("End of Trilinos details\n")

SET (CMAKE_CXX_FLAGS
     "-std=c++0x -Wall -DUSEH5FEDV2 -DIPPL_MPI -DMPICH_SKIP_MPICXX -DIPPL_DEBUG -DIPPL_DONT_POOL -DIPPL_USE_XDIV_RNG -DIPPL_LINUX -DIPPL_NO_STRINGSTREAM -DPETE_BITWISE_COPY -DIPPL_HAS_TEMPLATED_COMPLEX -DIPPL_USE_STANDARD_HEADERS -DIPPL_USE_PARTIAL_SPECIALIZATION -DIPPL_STDSTL -DIPPL_LONGLONG -DPARALLEL_IO -DNOCTAssert -Drestrict=__restrict__ -DH5_USE_16_API ${CMAKE_CXX_FLAGS}"
    )

MESSAGE (${CMAKE_CXX_FLAGS})

include_directories(
    ${CMAKE_CURRENT_SOURCE_PATH}
    ${CMAKE_CURRENT_SOURCE_PATH}/Commands
    ${GSL_INCLUDE_DIR}
    ${H5Hut_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
    ${TCLAP_INCLUDE_DIR}
    ${BOOST_INCLUDE_DIR}
    ${Trilinos_INCLUDE_DIRS}
    ${Trilinos_TPL_INCLUDE_DIRS}
    ${IPPL_INCLUDE_DIR}
)

link_directories(
   ${Trilinos_LIBRARY_DIRS}
   ${Trilinos_TPL_LIBRARY_DIRS}
   ${IPPL_LIBRARY_DIR}
)

SET (MaxwADIs_LIBS
    m
    ${GSL_LIBRARY}
    ${H5Hut_LIBRARY_C}
    ${H5Hut_LIBRARY}
    ${HDF5_LIBRARY}
    ${IPPL_LIBRARY}
    ${BOOST_REGEX_LIBRARY}
    ${BOOST_MPI_LIBRARY}
    ${BOOST_SERIALIZATION_LIBRARY}
    ${BOOST_TIMER_LIBRARY}
    ${BOOST_SYSTEM_LIBRARY}
    ${Metis_LIBRARY}
#    ${SZIP_LIBRARY}
    ${OPENSSL_LIBRARY}
    gfortran
    z
)

SET(CMAKE_BUILD_TYPE Debug)

add_subdirectory (test)
add_subdirectory (src)
