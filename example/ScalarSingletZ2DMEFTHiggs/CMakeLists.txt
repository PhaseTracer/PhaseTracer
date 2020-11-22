if(BUILD_WITH_FS_ScalarSingletZ2DMEFTHiggs)
  set(FS_version "v2.5.0")
  set(FS_model_name "ScalarSingletZ2DMEFTHiggs")
  set(SARAH_model_name "ScalarSingletZ2DM")

  find_package(FlexibleSUSY REQUIRED)
  find_package(GSL REQUIRED)

  add_executable(run_${FS_model_name} run_${FS_model_name}.cpp)

  target_include_directories(run_${FS_model_name} PRIVATE
    ${FlexibleSUSY_${FS_model_name}_INCLUDES}
    ${FlexibleSUSY_SM_INCLUDES}
    ${FlexibleSUSY_slha_INCLUDES}
    ${FlexibleSUSY_config_INCLUDES}
    ${FlexibleSUSY_src_INCLUDES}
  )

  target_link_libraries(run_${FS_model_name}
    ${FlexibleSUSY_${FS_model_name}}
    ${FlexibleSUSY_SM}
    ${FlexibleSUSY_src}
    GSL::gsl
  )

  set_target_properties(run_${FS_model_name}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
  )

endif(BUILD_WITH_FS_ScalarSingletZ2DMEFTHiggs)
