if(BUILD_WITH_FS)
  set(FS_version "v2.4.2")
  set(FS_model_name "THDMIISNMSSMBCsimple")
  set(SARAH_model_name "THDMS")

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

endif(BUILD_WITH_FS)
