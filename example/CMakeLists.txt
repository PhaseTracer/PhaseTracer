set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

link_libraries(
  phasetracer
  effectivepotential
)

add_executable(run_1D_test_model run_1D_test_model.cpp)
add_executable(run_2D_test_model run_2D_test_model.cpp)
add_executable(scan_Z2_scalar_singlet_model scan_Z2_scalar_singlet_model.cpp)
add_executable(time time.cpp)
add_executable(h_bar_expansion h_bar_expansion.cpp)
add_executable(run_xSM_msbar_noZ2_model run_xSM_msbar_noZ2_model.cpp)
add_executable(test_tree_ewsb test_tree_ewsb.cpp)
add_executable(test_file_reading test_file_reading.cpp)
add_executable(test_symmetry test_symmetry.cpp)
add_executable(run_SubcriticalTransitionExample run_SubcriticalTransitionExample.cpp)

# add any subdirs that might contain programs

macro(subdirlist result curdir)
  file(GLOB children RELATIVE ${curdir} ${curdir}/*)
  set(dirlist "")
  foreach(child ${children})
    if(IS_DIRECTORY ${curdir}/${child} AND EXISTS ${curdir}/${child}/CMakeLists.txt)
      list(APPEND dirlist ${child})
    endif()
  endforeach()
  set(${result} ${dirlist})
endmacro()

subdirlist(subdirs ${CMAKE_CURRENT_SOURCE_DIR})
foreach(subdir ${subdirs})
  add_subdirectory(${subdir})
endforeach()

