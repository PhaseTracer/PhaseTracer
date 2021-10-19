if(BUILD_WITH_BSMPT)

  find_package(BSMPT 2.3.2 QUIET)
  include(FetchContent)
  set(BSMPT_FETCHED_TARGETS OFF)
  if(NOT BSMPT_FOUND)
    find_package(Threads REQUIRED)
    FetchContent_Declare(
      BSMPT
      GIT_REPOSITORY https://github.com/phbasler/BSMPT.git
      GIT_TAG v2.3.2 )
    FetchContent_GetProperties(BSMPT)
    if(NOT BSMPT_POPULATED)
      set(BSMPT_FETCHED_TARGETS ON)
      option(UseLibCMAES OFF)
      FetchContent_Populate(BSMPT)
      add_subdirectory(${bsmpt_SOURCE_DIR} ${bsmpt_BINARY_DIR} )
    endif()
    set(BSMPT_FOUND TRUE)
  endif()

  add_executable(run_R2HDM run_R2HDM.cpp)

  add_executable(run_C2HDM run_C2HDM.cpp)

  add_executable(run_N2HDM run_N2HDM.cpp)

  if(BSMPT_FETCHED_TARGETS)
    target_link_libraries(run_R2HDM PRIVATE Minimizer Models Utility)
    target_link_libraries(run_C2HDM PRIVATE Minimizer Models Utility)
    target_link_libraries(run_N2HDM PRIVATE Minimizer Models Utility)
  else()
    target_link_libraries(run_R2HDM PRIVATE BSMPT::Minimizer BSMPT::Models BSMPT::Utility)
    target_link_libraries(run_C2HDM PRIVATE BSMPT::Minimizer BSMPT::Models BSMPT::Utility)
    target_link_libraries(run_N2HDM PRIVATE BSMPT::Minimizer BSMPT::Models BSMPT::Utility)
  endif()
  
endif(BUILD_WITH_BSMPT)
