# Find BubbleProfiler

set(BubbleProfiler "${PROJECT_SOURCE_DIR}/BubbleProfiler")

if(Git_FOUND)
  message("Git found: ${GIT_EXECUTABLE}")
elseif (NOT EXISTS ${BubbleProfiler} OR NOT EXISTS ${BubbleProfiler}/lib/libbubbler.a)
  message(FATAL_ERROR "Git not found.")
endif()

# Download BubbleProfiler if required
if(NOT EXISTS ${BubbleProfiler})
  message(STATUS "Downloading BubbleProfiler")
  set(BubbleProfiler_git_rep "https://github.com/bubbleprofiler/bubbleprofiler.git")
  execute_process(COMMAND git clone ${BubbleProfiler_git_rep} ${BubbleProfiler}
                 WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
	  )
endif()

# Build BubbleProfiler if required
if(NOT EXISTS ${BubbleProfiler}/lib/libbubbler.a)
  message(STATUS "Making BubbleProfiler")
  execute_process(COMMAND git checkout tags/v1.0.1
                 WORKING_DIRECTORY ${BubbleProfiler})
  execute_process(COMMAND cmake .
                  WORKING_DIRECTORY ${BubbleProfiler})
  execute_process(COMMAND ${CMAKE_MAKE_PROGRAM}
                 WORKING_DIRECTORY ${BubbleProfiler})
endif()

# find includes
find_path(BubbleProfiler_INCLUDES
  NAMES shooting.hpp
  PATHS ${BubbleProfiler}/include
)

# find libraries
find_library(BubbleProfiler_LIB
  NAMES libbubbler.a
  PATHS ${BubbleProfiler}/lib
)
