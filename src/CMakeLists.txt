file(GLOB SRC_FILES *.cpp)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
link_libraries(
  ${NLopt_LIBRARIES}
  ${Boost_LIBRARIES}
  ${EIGEN3_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
  ${ALGLIB_LIBRARIES}
  effectivepotential
)

if(BUILD_WITH_BP)
  link_libraries(
    ${BubbleProfiler_LIB}
    GSL::gsl
    GSL::gslcblas
  )
endif(BUILD_WITH_BP)

if(THREADS_HAVE_PTHREAD_ARG)
  target_compile_options(phasetracer "-pthread")
endif()

add_library(phasetracer SHARED ${SRC_FILES})
