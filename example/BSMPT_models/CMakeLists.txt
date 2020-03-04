if(BUILD_WITH_BSMPT)

  find_package(BSMPT REQUIRED)
  find_package(GSL REQUIRED)
  find_package(CMAES REQUIRED)
  find_package(OpenMP)

  link_libraries(
    ${BSMPT_MODELS_LIB}
    ${BSMPT_MINIMIZER_LIB}
    ${BSMPT_MODELS_LIB}
    ${CMAES_LIBRARY}
    GSL::gsl
    OpenMP::OpenMP_CXX
  )

  include_directories(
    ${BSMPT_MODELS_INCLUDES}
    ${BSMPT_MINIMIZER_INCLUDES}
    ${CMAES_INCLUDE_DIR}
  )

  add_executable(run_R2HDM run_R2HDM.cpp)
  add_executable(run_C2HDM run_C2HDM.cpp)
  add_executable(run_N2HDM run_N2HDM.cpp)
  
endif(BUILD_WITH_BSMPT)
