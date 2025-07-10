if(BUILD_WITH_DP)

find_package(OpenMP REQUIRED)

link_libraries(
   ${DeepPhase_LIB}
   OpenMP::OpenMP_CXX
)

add_executable(run_DP_xSM run_DP_xSM.cpp)

endif(BUILD_WITH_DP)