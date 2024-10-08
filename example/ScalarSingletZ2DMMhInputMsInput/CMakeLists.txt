if(BUILD_WITH_FS_ScalarSingletZ2DMMhInputMsInput)
  set(FS_version "v2.6.1")
  set(FS_model_name "ScalarSingletZ2DMMhInputMsInput")
  set(SARAH_model_name "ScalarSingletZ2DM")

  find_package(FlexibleSUSY REQUIRED)
  find_package(GSL REQUIRED)

  # mixed FS and by hand model so we get RGEs from FS + masses with singlet VEV
  add_executable(run_${FS_model_name}_withSingletVEVinPT run_${FS_model_name}_withSingletVEVinPT.cpp)

  target_include_directories(run_${FS_model_name}_withSingletVEVinPT PRIVATE
    ${FlexibleSUSY_${FS_model_name}_INCLUDES}
    ${FlexibleSUSY_SM_INCLUDES}
    ${FlexibleSUSY_slha_INCLUDES}
    ${FlexibleSUSY_config_INCLUDES}
    ${FlexibleSUSY_src_INCLUDES}
  )

  target_link_libraries(run_${FS_model_name}_withSingletVEVinPT
    ${FlexibleSUSY_${FS_model_name}}
    ${FlexibleSUSY_SM}
    ${FlexibleSUSY_src}
    GSL::gsl
    )
  
  set_target_properties(run_${FS_model_name}_withSingletVEVinPT
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
    )

  # mixed FS and by hand model so we get RGEs from FS + masses with singlet VEV
  add_executable(run_${FS_model_name}_withSingletVEVinPT_PRM run_${FS_model_name}_withSingletVEVinPT_PRM.cpp)

  target_include_directories(run_${FS_model_name}_withSingletVEVinPT_PRM PRIVATE
    ${FlexibleSUSY_${FS_model_name}_INCLUDES}
    ${FlexibleSUSY_SM_INCLUDES}
    ${FlexibleSUSY_slha_INCLUDES}
    ${FlexibleSUSY_config_INCLUDES}
    ${FlexibleSUSY_src_INCLUDES}
  )

  target_link_libraries(run_${FS_model_name}_withSingletVEVinPT_PRM
    ${FlexibleSUSY_${FS_model_name}}
    ${FlexibleSUSY_SM}
    ${FlexibleSUSY_src}
    GSL::gsl
    )
  
  set_target_properties(run_${FS_model_name}_withSingletVEVinPT_PRM
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
    )

endif(BUILD_WITH_FS_ScalarSingletZ2DMMhInputMsInput)
