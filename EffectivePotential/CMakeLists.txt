cmake_minimum_required(VERSION 3.1...3.15)
project(EffectivePotential VERSION 1.0.0 LANGUAGES CXX)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
enable_language(CXX)

# option for building example models
option(BUILD_WITH_EXAMPLE "Build example programs" OFF)

# packages

find_package(Eigen3 3.1 REQUIRED)
find_package(ALGLIB REQUIRED)

find_package(Boost 1.53.0 REQUIRED
  COMPONENTS
  filesystem
  system
  log
)
set(Boost_USE_MULTITHREADED ON)

# flags

set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# includes

include_directories(
  ${Boost_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  ${ALGLIB_INCLUDES} 
  ${PROJECT_SOURCE_DIR}/include/effectivepotential
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)

# add relevant dirs

add_subdirectory(src)
add_subdirectory(example)
