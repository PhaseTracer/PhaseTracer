cmake_minimum_required(VERSION 3.9)
project(PhaseTracer VERSION 2.1.0 LANGUAGES CXX)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
enable_language(CXX)
enable_testing()

# option for building models with FlexibleSUSY
option(BUILD_WITH_FS "Build FlexibleSUSY example program (THDMS with NMSSM BCs) " OFF)

option(BUILD_WITH_FS_ScalarSingletZ2DM "Build FlexibleSUSY ScalarSingletZ2DM program" OFF)

option(BUILD_WITH_FS_ScalarSingletZ2DMEWSBoutputlamHEFTHiggs "Build FlexibleSUSY_ScalarSingletZ2DMEWSBoutputlamHEFTHiggs program" OFF)

# option for building models with BSMPT
option(BUILD_WITH_BSMPT "Build BSMPT example programs" OFF)

# option for including BubbleProfiler
option(BUILD_WITH_BP "Build example programs including BubbleProfiler" OFF)

# packages

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 OLD)
endif()

find_package(Boost 1.53.0 REQUIRED
  COMPONENTS
  filesystem
  system
  log
)
set(Boost_USE_MULTITHREADED ON)

find_package(ALGLIB REQUIRED)
find_package(Eigen3 3.1 REQUIRED)
find_package(NLopt REQUIRED)
find_package(Threads REQUIRED)
find_package(Git)
# find_package(EffPotential REQUIRED)
add_subdirectory(EffectivePotential)

if(BUILD_WITH_BP)
  add_definitions(-DBUILD_WITH_BP)
  find_package(GSL REQUIRED)
  find_package(BubbleProfiler REQUIRED)
endif(BUILD_WITH_BP)

# flags

set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# includes

get_target_property(EFFPOTENTIAL_INCLUDE_DIRS effectivepotential INTERFACE_INCLUDE_DIRECTORIES)  
include_directories(
  ${Boost_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  ${NLopt_INCLUDE_DIRS}
  ${EFFPOTENTIAL_INCLUDE_DIRS}
  ${ALGLIB_INCLUDES}
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)

if(BUILD_WITH_BP)
  include_directories(
    ${PROJECT_SOURCE_DIR}
  )
endif(BUILD_WITH_BP)

# add relevant dirs

add_subdirectory(src)
add_subdirectory(example)
add_subdirectory(unit_tests)
