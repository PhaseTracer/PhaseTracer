# ============================================================================
#  Standalone build for the `pt_bounce` Python module (BubbleDet validation).
#
#  Links against the already-built shared libraries
#      <PT_ROOT>/EffectivePotential/lib/libeffectivepotential.so
#      <PT_ROOT>/lib/libphasetracer.so
#  (both built SHARED), so build the main project first.
#
#  Build:
#      cd example/BubbleDetBridge
#      cmake -S . -B build
#      cmake --build build -j
#  Produces  build/pt_bounce.<abi>.so  which Python can import.
# ============================================================================
cmake_minimum_required(VERSION 3.18)
project(pt_bounce_bridge LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# This dir is <PT_ROOT>/example/BubbleDetBridge
get_filename_component(PT_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
set(PT_ROOT "${PT_ROOT}" CACHE PATH "PhaseTracer source root")
list(APPEND CMAKE_MODULE_PATH "${PT_ROOT}/cmake")

message(STATUS "PhaseTracer root: ${PT_ROOT}")

# --- dependencies ----------------------------------------------------------
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(Boost 1.53.0 REQUIRED)
find_package(Eigen3 3.1 REQUIRED)
find_package(ALGLIB REQUIRED)
find_package(GSL REQUIRED)
find_package(NLopt REQUIRED)

find_package(pybind11 CONFIG QUIET)
if(NOT pybind11_FOUND)
  message(STATUS "pybind11 not found - fetching it via FetchContent")
  include(FetchContent)
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v2.13.6
  )
  FetchContent_MakeAvailable(pybind11)
else()
  message(STATUS "Using installed pybind11: ${pybind11_DIR}")
endif()

# --- the module ------------------------------------------------------------
pybind11_add_module(pt_bounce bounce_bridge.cpp)

target_include_directories(pt_bounce PRIVATE
  "${PT_ROOT}/EffectivePotential/include"                      # models/..., effectivepotential/...
  "${PT_ROOT}/EffectivePotential/include/effectivepotential"   # bare potential.hpp, pow.hpp
  "${PT_ROOT}/include"                                         # phasetracer.hpp, property.hpp
  "${PT_ROOT}"
  ${Boost_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  ${ALGLIB_INCLUDES}
  ${NLopt_INCLUDE_DIRS}
)

target_link_libraries(pt_bounce PRIVATE
  "${PT_ROOT}/EffectivePotential/lib/libeffectivepotential.so"
  "${PT_ROOT}/lib/libphasetracer.so"
  ${ALGLIB_LIBRARIES}
  GSL::gsl
  GSL::gslcblas
  ${NLopt_LIBRARIES}
)

# Make the built libs findable at import time without setting LD_LIBRARY_PATH.
set_target_properties(pt_bounce PROPERTIES
  BUILD_RPATH "${PT_ROOT}/EffectivePotential/lib;${PT_ROOT}/lib"
)
