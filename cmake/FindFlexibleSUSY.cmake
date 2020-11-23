# Find FlexibleSUSY and a particular model
# The model is defined by FS_model_name and SARAH_model_name

set(FS "${PROJECT_SOURCE_DIR}/FlexibleSUSY")

# Download FlexibleSUSY if required
if(NOT EXISTS ${FS})
  message(STATUS "Downloading FlexibleSUSY")
  set(FS_git_rep "https://github.com/FlexibleSUSY/FlexibleSUSY.git")
  execute_process(COMMAND git clone ${FS_git_rep} ${FS}
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
		  )
endif()

if(NOT DEFINED FS_PARALLEL_MAKE)
  set(FS_PARALLEL_MAKE "1")
endif()

# Compile model
if(NOT EXISTS ${FS}/models/${FS_model_name}/lib${FS_model_name}.a)
  if(NOT EXISTS ${FS}/models/${FS_model_name})
    message(STATUS "Creating ${FS_model_name} model files")
    execute_process(COMMAND git checkout tags/${FS_version}
                    WORKING_DIRECTORY ${FS}
                   )
    execute_process(COMMAND cp -r ${FS_model_name} ${FS}/models/
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/example/${FS_model_name}/FS_generated_code/
    )
  endif()

  message(STATUS "Making ${FS_model_name} model")
  execute_process(COMMAND ./configure --with-models=${FS_model_name} --disable-threads --disable-meta
                  WORKING_DIRECTORY ${FS}
  )
  	
  execute_process(COMMAND ${CMAKE_MAKE_PROGRAM} -j${FS_PARALLEL_MAKE}
                    WORKING_DIRECTORY ${FS}
  )
 
endif()

# find includes
find_path(FlexibleSUSY_config_INCLUDES
  NAMES config.h
  PATHS ${FS}/config
)

find_path(FlexibleSUSY_SM_INCLUDES
  NAMES standard_model.hpp
  PATHS ${FS}/model_specific/SM
)

find_path(FlexibleSUSY_src_INCLUDES
  NAMES problems.hpp
  PATHS ${FS}/src
)

find_path(FlexibleSUSY_slha_INCLUDES
  NAMES slhaea.h
  PATHS ${FS}/slhaea
)

find_path(FlexibleSUSY_${FS_model_name}_INCLUDES
  NAMES ${FS_model_name}_input_parameters.hpp
  PATHS 
  ${FS}/models/${FS_model_name}
)

# find libraries
find_library(FlexibleSUSY_src
  NAMES libflexisusy.a
  PATHS ${FS}/src
)

find_library(FlexibleSUSY_SM
  NAMES libmodel_specific_SM.a 
  PATHS ${FS}/model_specific/SM
)

find_library(FlexibleSUSY_${FS_model_name}
  NAMES lib${FS_model_name}.a
  PATHS ${FS}/models/${FS_model_name}
)
