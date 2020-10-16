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
  execute_process(COMMAND git checkout tags/v2.0.1
                  WORKING_DIRECTORY ${BSMPT})
  execute_process(COMMAND cmake -DUseLibCMAES=OFF .
                  WORKING_DIRECTORY ${BSMPT})
  execute_process(COMMAND ${CMAKE_MAKE_PROGRAM}
                  WORKING_DIRECTORY ${BSMPT})
endif()

# find includes
find_path(BSMPT_MODELS_INCLUDES
  NAMES ClassPotentialOrigin.h
  PATHS ${BSMPT}/include/BSMPT/models
)

find_path(BSMPT_MINIMIZER_INCLUDES
  NAMES Minimizer.h
  PATHS ${BSMPT}/include/BSMPT/minimizer
)

find_path(BSMPT_CONFIG_INCLUDES
  NAMES BSMPT/config.h 
  PATHS ${BSMPT}/include
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

find_library(BSMPT_THERMAL_LIB
  NAMES libThermalFunctions.a
  PATHS ${BSMPT}/src/ThermalFunctions
)

find_library(BSMPT_UTILITY_LIB
  NAMES libUtility.a
  PATHS ${BSMPT}/src
)

