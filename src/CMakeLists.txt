file(GLOB SRC_FILES *.cpp)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
link_libraries(
  ${NLopt_LIBRARIES}
  ${Boost_LIBRARIES}
  ${EIGEN3_LIBRARIES}
  effectivepotential
)
add_library(phasetracer SHARED ${SRC_FILES})
