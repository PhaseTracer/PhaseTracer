add_executable(run_false_vacuum_decay_rate run_false_vacuum_decay_rate.cpp)
add_executable(run_equation_of_state run_equation_of_state.cpp)
add_executable(run_friedmann_evolution run_friedmann_evolution.cpp)
add_executable(run_thermo_finder run_thermo_finder.cpp)

# To benchmark the action calculator speed up we call the OpenMP api, so it 
# needs to be linked here.
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(run_false_vacuum_decay_rate PRIVATE OpenMP::OpenMP_CXX)
endif()

