file(GLOB SRC_FILES *.cpp)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

if(THREADS_HAVE_PTHREAD_ARG)
  target_compile_options(phasetracer "-pthread")
endif()

add_library(phasetracer SHARED ${SRC_FILES})

target_link_libraries(phasetracer PRIVATE
  ${NLopt_LIBRARIES}
  ${Boost_LIBRARIES}
  ${EIGEN3_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
  ${ALGLIB_LIBRARIES}
  effectivepotential
)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(phasetracer PUBLIC OpenMP::OpenMP_CXX)
else()
  message(WARNING
    "OpenMP not found: the bounce-action loop will run serially.")
endif()

if(BUILD_WITH_BP)
  target_link_libraries(phasetracer PRIVATE
    ${BubbleProfiler_LIB}
    GSL::gsl
    GSL::gslcblas
  )
endif(BUILD_WITH_BP)

if(BUILD_WITH_HG)
  target_link_libraries(phasetracer PRIVATE
    ${HydroGrav_LIBRARIES}
    GSL::gsl
    GSL::gslcblas
  )
endif(BUILD_WITH_HG)
