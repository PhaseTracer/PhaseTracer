set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

link_libraries(
  phasetracer
  effectivepotential
)

add_executable(run_1D_test_model run_1D_test_model.cpp)
add_executable(run_2D_test_model run_2D_test_model.cpp)
add_executable(scan_Z2_scalar_singlet_model scan_Z2_scalar_singlet_model.cpp)
add_executable(time time.cpp)

add_executable(run_ToyModel run_ToyModel.cpp)
add_executable(run_RSS run_RSS.cpp)
add_executable(run_supercool run_supercool.cpp)
add_executable(run_ConcurrentThreePhaseTransition_scaled run_ConcurrentThreePhaseTransition_scaled.cpp)
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

