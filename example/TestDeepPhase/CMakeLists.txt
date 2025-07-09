if(BUILD_WITH_DP)
    link_libraries(
        ${DeepPhase_LIB}
    )
    add_executable(run_DP_xSM run_DP_xSM.cpp)
endif(BUILD_WITH_DP)