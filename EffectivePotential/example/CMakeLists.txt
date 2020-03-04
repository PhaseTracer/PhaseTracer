if(BUILD_WITH_EXAMPLE)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

link_libraries(
  effectivepotential
  ${EIGEN3_LIBRARIES}
  ${ALGLIB_LIBRARIES} 
)

add_executable(simplified_2d simplified_2d.cpp)
add_executable(Z2_scalar_singlet Z2_scalar_singlet.cpp)
endif(BUILD_WITH_EXAMPLE)

