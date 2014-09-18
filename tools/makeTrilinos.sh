CC=mpicc
CXX=mpicxx
CPP="mpicxx -E"
F77=mpif77

# METIS
METIS_PATH=/path/to/metis-5.0.2/build
METIS_INCLUDE_PATH=$METIS_PATH/include
METIS_LIBRARY_PATH=$METIS_PATH/lib

# PARMETIS
PARMETIS_PATH=/path/to/parmetis-4.0.2/build
PARMETIS_INCLUDE_PATH=$PARMETIS_PATH/include
PARMETIS_LIBRARY_PATH=$PARMETIS_PATH/lib

SCALAPACK_INCLUDE_PATH=/path/to/scalapack/

SUPERLU_PATH=/path/to/SuperLU_DIST_3.0/
SUPERLU_INCLUDE_PATH=$SUPERLU_PATH/SRC
SUPERLU_LIBRARY_PATH=$SUPERLU_PATH/lib

TRILINOS_PATH=/path/to/trilinos-10.8.3/src
cmake \
           --prefix=${TRILINOS_PATH}/build \
           -DCMAKE_INSTALL_PREFIX:PATH=${TRILINOS_PATH}/build \
           -DCMAKE_CXX_FLAGS:STRING="-DMPICH_IGNORE_CXX_SEEK -fPIC" \
           -DCMAKE_C_FLAGS:STRING="-DMPICH_IGNORE_CXX_SEEK -fPIC" \
           -DCMAKE_Fortran_FLAGS:STRING="-fPIC" \
           -D CMAKE_BUILD_TYPE:STRING=DEBUG \
           -D TPL_ENABLE_SuperLUDist:BOOL=ON\
           -D TPL_ENABLE_MPI:BOOL=ON \
           -D TPL_ENABLE_BLAS:BOOL=ON \
           -D TPL_ENABLE_LAPACK:BOOL=ON \
           -D TPL_ENABLE_ParMETIS:BOOL=ON \
           -D TPL_ENABLE_METIS:BOOL=OFF \
           -D ParMETIS_INCLUDE_DIRS:PATH="${PARMETIS_INCLUDE_PATH}" \
           -D ParMETIS_LIBRARY_DIRS:PATH="${PARMETIS_LIBRARY_PATH}" \
           -D ParMETIS_LIBRARY_NAMES:STRING="parmetis" \
           -D ParMETIS_LIBRARIES:STRING="libparmetis.a" \
	   -D METIS_INCLUDE_DIRS:PATH="${METIS_INCLUDE_PATH}" \
	   -D METIS_LIBRARY_DIRS:PATH="${METIS_LIBRARY_PATH}" \
	   -D METIS_LIBRARY_NAMES:STRING="metis" \
	   -D METIS_LIBRARIES:STRING="metis" \
           -D SuperLUDist_INCLUDE_DIRS:FILEPATH="${SUPERLU_INCLUDE_PATH}" \
           -D SuperLUDist_LIBRARY_DIRS:FILEPATH="${SUPERLU_LIBRARY_PATH}" \
           -D SuperLUDist_LIBRARY_NAMES:STRING="superlu_dist_3.0" \
           -D SuperLUDist_LIBRARIES="superlu_dist_3.0" \
           -D TPL_SuperLUDist_LIBRARIES="${SUPERLU_LIBRARY_PATH}/libsuperlu_dist_3.0.a" \
           -D SCALAPACK_INCLUDE_DIRS:FILEPATH="${SCALAPACK_INCLUDE_PATH}" \
           -D SCALAPACK_LIBRARY_DIRS:FILEPATH="${SCALAPACK_INCLUDE_PATH}" \
           -D SCALAPACK_LIBRARY_NAMES:STRING="scalapack" \
           -D Trilinos_ENABLE_Belos:BOOL=ON \
           -D Trilinos_ENABLE_Epetra:BOOL=ON \
           -D Trilinos_ENABLE_EpetraExt:BOOL=ON \
           -D Trilinos_ENABLE_Ifpack:BOOL=ON \
           -D Trilinos_ENABLE_ML:BOOL=ON \
           -D Trilinos_ENABLE_Amesos:BOOL=ON \
           -D Amesos_ENABLE_BLACS:BOOL=ON \
           -D Amesos_ENABLE_SuperLUDist:BOOL=ON \
           -D Amesos_ENABLE_SCALAPACK:BOOL=ON \
           -D Trilinos_ENABLE_Amesos-superlu:BOOL=ON\
           -D Trilinos_ENABLE_AztecOO:BOOL=ON \
           -D Trilinos_ENABLE_Teuchos:BOOL=ON \
           -D Trilinos_ENABLE_Aztecoo-Teuchos:BOOL=ON \
           -D Trilinos_ENABLE_Teuchos-Extended:BOOL=ON \
           -D Trilinos_ENABLE_Isorropia:BOOL=ON \
           -D Trilinos_ENABLE_Isorropia-Epetraext:BOOL=ON \
           -D Trilinos_ENABLE_Didasko:BOOL=OFF \
           -D Didasko_ENABLE_TESTS=OFF \
           -D Didasko_ENABLE_EXAMPLES=OFF \
           -D Trilinos_ENABLE_TESTS:BOOL=OFF \
           $EXTRA_ARGS \
           ${TRILINOS_PATH}
