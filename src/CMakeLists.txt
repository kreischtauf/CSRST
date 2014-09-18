configure_file(${CMAKE_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/csrst_config.h)

SET (MaxwADIs_SRCS
    TETMSolver.cpp
    Bend.cpp
    Distribution.cpp
    PartBunch.cpp
    Physics.cpp
    utils.cpp
    PoissonSolver.cpp
    ZerothOrderShapeFunction.cpp
    FirstOrderShapeFunction.cpp
    SecondOrderShapeFunction.cpp
    Communicator.cpp
    Commands/AllocateFields.cpp
    Commands/ASCIIVTKFile.cpp
    Commands/BinaryVTKFile.cpp
    Commands/CalcEnergy.cpp
    Commands/LoadH5Part.cpp
    Commands/PlotASCIIVTK.cpp
    Commands/PlotBinaryVTK.cpp
    Commands/SaveH5Part.cpp
    Commands/DataSink.cpp
    Commands/SaveEMFields.cpp
    Commands/IntegrateElong.cpp
    gauss_legendre.c
)

add_library( maxwadis ${MaxwADIs_SRCS} )

add_executable (csrst main.cpp)
target_link_libraries( csrst maxwadis ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${MaxwADIs_LIBS})
