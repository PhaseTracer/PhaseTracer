cmake_minimum_required(VERSION 3.9)
project(PhaseTracer VERSION 1.0.0 LANGUAGES CXX)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
enable_language(CXX)

# option for building models with FlexibleSUSY
option(BUILD_WITH_FS "Build FlexibleSUSY example programs" OFF)

# option for building models with BSMPT
option(BUILD_WITH_BSMPT "Build BSMPT example programs" OFF)

# packages

find_package(Boost 1.53.0 REQUIRED
  COMPONENTS
  filesystem
  log
)
find_package(Eigen3 3.1 REQUIRED)
find_package(NLopt REQUIRED)
find_package(Threads REQUIRED)
# find_package(EffPotential REQUIRED)
add_subdirectory(EffectivePotential)

# flags

set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# includes

get_target_property(EFFPOTENTIAL_INCLUDE_DIRS effectivepotential INTERFACE_INCLUDE_DIRECTORIES)  
include_directories(
  ${Boost_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  ${NLopt_INCLUDE_DIRS}
  ${EFFPOTENTIAL_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)

# add relevant dirs

add_subdirectory(src)
add_subdirectory(example)
add_subdirectory(unit_tests)
