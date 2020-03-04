# Find BSMPT

set(BSMPT "${PROJECT_SOURCE_DIR}/BSMPT")

# Download BSMPT if required
if(NOT EXISTS ${BSMPT})
  message(STATUS "Downloading BSMPT")
  set(BSMPT_git_rep "https://github.com/phbasler/BSMPT.git")
  execute_process(COMMAND git clone ${BSMPT_git_rep} ${BSMPT}
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
		  )
endif()

# Build BSMPT if required
if(NOT EXISTS ${BSMPT}/src/models/libModels.a)
  message(STATUS "Making BSMPT")
  execute_process(COMMAND cmake .
                  WORKING_DIRECTORY ${BSMPT})
  execute_process(COMMAND ${CMAKE_MAKE_PROGRAM}
                  WORKING_DIRECTORY ${BSMPT})
endif()

# find includes
find_path(BSMPT_MODELS_INCLUDES
  NAMES ClassPotentialOrigin.h
  PATHS ${BSMPT}/src/models
)

find_path(BSMPT_MINIMIZER_INCLUDES
  NAMES Minimizer.h
  PATHS ${BSMPT}/src/minimizer
)


# find libraries
find_library(BSMPT_MODELS_LIB
  NAMES libModels.a
  PATHS ${BSMPT}/src/models
)

find_library(BSMPT_MINIMIZER_LIB
  NAMES libMinimizer.a
  PATHS ${BSMPT}/src/minimizer
)
