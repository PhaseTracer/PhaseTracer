file(GLOB SRC_FILES *.cpp)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
link_libraries(
  ${Boost_LIBRARIES}
  ${EIGEN3_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
  effectivepotential
)

if(TARGET NLopt::nlopt)
  link_libraries(NLopt::nlopt)
else()
  link_libraries(NLopt::nlopt_cxx)
endif()

if(THREADS_HAVE_PTHREAD_ARG)
  target_compile_options(phasetracer "-pthread")
endif()

add_library(phasetracer SHARED ${SRC_FILES})
