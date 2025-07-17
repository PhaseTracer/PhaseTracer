# Find DeepPhase

set(DeepPhase_GIT_BRANCH "main" CACHE STRING "Branch of DeepPhase to clone")
set(DeepPhase "${PROJECT_SOURCE_DIR}/DeepPhase")

if(Git_FOUND)
  message("Git found: ${GIT_EXECUTABLE}")
elseif (NOT EXISTS ${DeepPhase} OR NOT EXISTS ${DeepPhase}/lib/libDeepPhaseLib.a)
  message(FATAL_ERROR "Git not found.")
endif()

# Download DeepPhase if required
if(NOT EXISTS ${DeepPhase})
  message(STATUS "Downloading DeepPhase (branch: ${DeepPhase_GIT_BRANCH})")
  set(DeepPhase_git_rep "https://github.com/Comparison-Project/DeepPhase.git")
  execute_process(COMMAND git clone --branch ${DeepPhase_GIT_BRANCH} --single-branch ${DeepPhase_git_rep} ${DeepPhase}
                 WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
endif()

# Build DeepPhase if required
if(NOT EXISTS ${DeepPhase}/lib/libDeepPhaseLib.a)
  message(STATUS "Making DeepPhase")
  # execute_process(COMMAND git checkout tags/v1.0.0
  #                WORKING_DIRECTORY ${DeepPhase})
  execute_process(COMMAND cmake -DCMAKE_CXX_FLAGS="-fPIC" -DBUILD_WITH_UNIT_TESTS=OFF -DBUILD_WITH_EXAMPLES=OFF .
                  WORKING_DIRECTORY ${DeepPhase})
  execute_process(COMMAND ${CMAKE_MAKE_PROGRAM}
                 WORKING_DIRECTORY ${DeepPhase})
endif()

# find includes
find_path(DeepPhase_INCLUDES
  NAMES deepphase.hpp
  PATHS ${DeepPhase}/include
)

# find libraries
find_library(DeepPhase_LIB
  NAMES libDeepPhaseLib.a
  PATHS ${DeepPhase}/lib
)