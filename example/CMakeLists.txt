set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

link_libraries(
  phasetracer
  effectivepotential
)

add_executable(run_1D_test_model run_1D_test_model.cpp)
add_executable(run_2D_test_model run_2D_test_model.cpp)
add_executable(scan_Z2_scalar_singlet_model scan_Z2_scalar_singlet_model.cpp)
add_executable(time time.cpp)
if(TARGET NLopt::nlopt)
  target_link_libraries(run_1D_test_model PUBLIC NLopt::nlopt)
  target_link_libraries(run_2D_test_model PUBLIC NLopt::nlopt)
  target_link_libraries(scan_Z2_scalar_singlet_model PUBLIC NLopt::nlopt)
  target_link_libraries(time PUBLIC NLopt::nlopt)
else()
  target_link_libraries(run_1D_test_model PUBLIC NLopt::nlopt_cxx)
  target_link_libraries(run_2D_test_model PUBLIC NLopt::nlopt_cxx)
  target_link_libraries(scan_Z2_scalar_singlet_model PUBLIC NLopt::nlopt_cxx)
  target_link_libraries(time PUBLIC NLopt::nlopt_cxx)
endif()

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

