message("adding test SCATTERGATHER found in ${CMAKE_CURRENT_SOURCE_DIR}")

include_directories (
    ${CMAKE_SOURCE_DIR}/src
)

link_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}
)

add_executable (ScatterGather EXCLUDE_FROM_ALL ../src/Physics.cpp ../src/utils.cpp scattergather.cpp)
target_link_libraries (ScatterGather ippl ${MaxwADIs_LIBS})

add_executable (calcCurrent EXCLUDE_FROM_ALL ../src/utils.cpp calcCurrent.cpp)
target_link_libraries (calcCurrent ippl ${MaxwADIs_LIBS})

add_executable (calcCurrentS1 EXCLUDE_FROM_ALL ../src/utils.cpp calcCurrentS1.cpp)
target_link_libraries (calcCurrentS1 maxwadis ippl ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${MaxwADIs_LIBS})
