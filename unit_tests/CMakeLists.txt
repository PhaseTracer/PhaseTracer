file(GLOB TEST_FILES test*.cpp)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

find_package(NLopt REQUIRED)

link_libraries(
  phasetracer
  effectivepotential
  ${NLopt_LIBRARIES}
)

add_executable(unit_tests catch_main.cpp ${TEST_FILES})

# add tests

add_test(NAME UnitTests COMMAND unit_tests)
