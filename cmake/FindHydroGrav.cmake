# Find HydroGrav

set(HydroGrav_GIT_BRANCH "main" CACHE STRING "Branch of HydroGrav to clone")
set(HydroGrav "${PROJECT_SOURCE_DIR}/HydroGrav")

if(Git_FOUND)
  message("Git found: ${GIT_EXECUTABLE}")
elseif (NOT EXISTS ${HydroGrav} OR NOT EXISTS ${HydroGrav}/lib/libHydroGravLib.a)
  message(FATAL_ERROR "Git not found.")
endif()

# Download HydroGrav if required
if(NOT EXISTS ${HydroGrav})
  message(STATUS "Downloading HydroGrav (branch: ${HydroGrav_GIT_BRANCH})")
  set(HydroGrav_git_rep "https://github.com/Hydro-Grav/HydroGrav.git")
  execute_process(COMMAND git clone --branch ${HydroGrav_GIT_BRANCH} --single-branch ${HydroGrav_git_rep} ${HydroGrav}
                 WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
endif()

# Build HydroGrav if required
if(NOT EXISTS ${HydroGrav}/lib/libHydroGravLib.a)
  message(STATUS "Making HydroGrav")
  # execute_process(COMMAND git checkout tags/v1.0.0
  #                WORKING_DIRECTORY ${HydroGrav})
  execute_process(COMMAND cmake -DCMAKE_CXX_FLAGS="-fPIC" -DBUILD_WITH_UNIT_TESTS=OFF -DBUILD_WITH_EXAMPLES=OFF .
                  WORKING_DIRECTORY ${HydroGrav})
  execute_process(COMMAND ${CMAKE_MAKE_PROGRAM}
                 WORKING_DIRECTORY ${HydroGrav})
endif()

# find includes
find_path(HydroGrav_INCLUDES
  NAMES hydrograv.hpp
  PATHS ${HydroGrav}/include
)

# find libraries
find_library(HydroGrav_LIB
  NAMES libHydroGravLib.a
  PATHS ${HydroGrav}/lib
)

# Find HydroGrav's dependencies (FINUFFT and FFTW3)
# These are required because HydroGrav is a static library
# find_package(GSL REQUIRED)
# find_library(FINUFFT_LIBRARY finufft PATHS /usr/local/lib REQUIRED)
# find_library(FFTW3_LIBRARY fftw3 REQUIRED)
# find_library(FFTW3F_LIBRARY fftw3f REQUIRED)
# find_library(FFTW3_OMP_LIBRARY fftw3_omp REQUIRED)
# find_library(FFTW3F_OMP_LIBRARY fftw3f_omp REQUIRED)

# Package all HydroGrav-related libraries
set(HydroGrav_LIBRARIES
  ${HydroGrav_LIB}
  GSL::gsl
  GSL::gslcblas
  # ${FINUFFT_LIBRARY}
  # ${FFTW3_LIBRARY}
  # ${FFTW3F_LIBRARY}
  # ${FFTW3_OMP_LIBRARY}
  # ${FFTW3F_OMP_LIBRARY}
)