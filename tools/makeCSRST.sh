#!/bin/bash

export CSRST_ROOT=/path/to/CSRST

export HDF5_ROOT=/path/to/hdf5
export HDF5_INCLUDE_PATH=$HDF5_ROOT/include
export HDF5_LIBRARY_PATH=$HDF5_ROOT/lib

export H5hut=/path/to/H5hut

export IPPL_ROOT=/path/to/ippl-edge
export IPPL_PREFIX=${IPPL_ROOT}/build

export Metis_ROOT=/path/to/metis-5.0.2/build
export TRILINOS_INCLUDE_PATH=/path/to/trilinos-10.8.3/build

export TCLAP_INCLUDE_PATH=/path/to/tclap-1.2.1/include
rm CMakeCache.txt

CXX=mpicxx cmake \
	-DCMAKE_CXX_FLAGS="-pg -fpermissive -Wno-literal-suffix -Wno-unused-local-typedefs" \
	-DCMAKE_EXE_LINKER_FLAGS="-pg" \
	$CSRST_ROOT
