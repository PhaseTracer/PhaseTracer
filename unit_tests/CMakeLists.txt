file(GLOB TEST_FILES test*.cpp)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

link_libraries(
  phasetracer
  effectivepotential
)

add_executable(unit_tests catch_main.cpp ${TEST_FILES})
if(TARGET NLopt::nlopt)
  target_link_libraries(unit_tests PUBLIC NLopt::nlopt)
else()
  target_link_libraries(unit_tests PUBLIC NLopt::nlopt_cxx)
endif()
