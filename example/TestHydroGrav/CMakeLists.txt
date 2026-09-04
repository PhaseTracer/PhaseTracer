if(BUILD_WITH_HG)

find_package(OpenMP REQUIRED)

add_executable(run_hydrograv run_hydrograv.cpp)

set(_hydrogravscripts_link_deps
	-Wl,--start-group
	${DeepPhase_LIBRARIES}
	${HydroGrav_LIBRARIES}
	${NLopt_LIBRARIES}
	OpenMP::OpenMP_CXX
	GSL::gsl
	GSL::gslcblas
	-Wl,--end-group
)

target_link_libraries(run_hydrograv PRIVATE ${_hydrogravscripts_link_deps})

endif(BUILD_WITH_HG)