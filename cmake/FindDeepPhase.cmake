# Find DeepPhase

set(DeepPhase "${PROJECT_SOURCE_DIR}/DeepPhase")

if(Git_FOUND)
  message("Git found: ${GIT_EXECUTABLE}")
elseif (NOT EXISTS ${DeepPhase} OR NOT EXISTS ${DeepPhase}/lib/libDeepPhaseLib.a)
  message(FATAL_ERROR "Git not found.")
endif()

# Download DeepPhase if required
if(NOT EXISTS ${DeepPhase})
  message(STATUS "Downloading DeepPhase")
  set(DeepPhase_git_rep "https://github.com/William-Searle/DeepPhase.git")
  execute_process(COMMAND git clone ${DeepPhase_git_rep} ${DeepPhase}
                 WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
	  )
endif()

# Build DeepPhase if required
if(NOT EXISTS ${DeepPhase}/lib/libDeepPhaseLib.a)
  message(STATUS "Making DeepPhase")
  # execute_process(COMMAND git checkout tags/v1.0.0
  #                WORKING_DIRECTORY ${DeepPhase})
  execute_process(COMMAND cmake -DCMAKE_CXX_FLAGS="-fPIC" .
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