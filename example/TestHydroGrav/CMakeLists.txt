# These examples use the HydroGrav bridge directly
if(BUILD_WITH_HG)
  add_executable(run_hydrograv run_hydrograv.cpp)
  add_executable(run_easy_hydrograv run_easy_hydrograv.cpp)
endif()
