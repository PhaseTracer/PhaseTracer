# Python interface: the `phasetracer` package, built around the extension module
# phasetracer._phasetracer.
#
#   Development build:  cmake -DBUILD_PYTHON=ON .. && make _phasetracer
#                       then PYTHONPATH=<build>/python python3 -c "import phasetracer"
#   pip:                pip install .   (from the repository root; see pyproject.toml)

# An extension module only needs Development.Module. A plain build asks for Development, the same
# as other subprojects may already have done: FindPython caches its results, and asking again for
# a different component set fails.
if(DEFINED SKBUILD)
  find_package(Python 3.8 COMPONENTS Interpreter Development.Module REQUIRED)
else()
  find_package(Python 3.8 COMPONENTS Interpreter Development REQUIRED)
endif()

# pybind11: one already loaded by another subproject, an installed copy (e.g. `pip install
# pybind11`, or pip's build environment), or a download
if(COMMAND pybind11_add_module)
  message(STATUS "Using the pybind11 already loaded by another subproject")
else()
  if(NOT pybind11_DIR)
    execute_process(
      COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
      OUTPUT_VARIABLE PT_PYBIND11_CMAKEDIR
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
    )
    if(PT_PYBIND11_CMAKEDIR)
      set(pybind11_DIR "${PT_PYBIND11_CMAKEDIR}")
    endif()
  endif()
  find_package(pybind11 CONFIG QUIET)
  if(pybind11_FOUND)
    message(STATUS "Using pybind11 ${pybind11_VERSION}: ${pybind11_DIR}")
  else()
    message(STATUS "pybind11 not found - fetching it")
    include(FetchContent)
    FetchContent_Declare(
      pybind11
      GIT_REPOSITORY https://github.com/pybind/pybind11.git
      GIT_TAG        v2.13.6
    )
    FetchContent_MakeAvailable(pybind11)
  endif()
endif()

file(GLOB PT_PYTHON_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
pybind11_add_module(_phasetracer ${PT_PYTHON_SOURCES})
target_link_libraries(_phasetracer PRIVATE phasetracer effectivepotential)

# A static Boost.Log (from scripts/install_dependencies.sh --local) lives inside
# libphasetracer.so; linking it again would give the module its own logging core.
if(NOT Boost_USE_STATIC_LIBS)
  target_link_libraries(_phasetracer PRIVATE ${Boost_LIBRARIES})
endif()

if(DEFINED SKBUILD)
  # pip: the wheel carries the module and the shared libraries it needs, side by side
  if(APPLE)
    set(PT_PYTHON_RPATH "@loader_path")
  else()
    set(PT_PYTHON_RPATH "$ORIGIN")
  endif()
  set_target_properties(_phasetracer phasetracer effectivepotential PROPERTIES INSTALL_RPATH "${PT_PYTHON_RPATH}")
  install(TARGETS _phasetracer phasetracer effectivepotential
          LIBRARY DESTINATION phasetracer COMPONENT python)

  # ALGLIB built by --local is a shared library that is not on the system search path
  get_filename_component(PT_ALGLIB_DIR "${ALGLIB_LIBRARIES}" DIRECTORY)
  if(EXISTS "${PT_DEPS_PREFIX}/lib" AND PT_ALGLIB_DIR STREQUAL "${PT_DEPS_PREFIX}/lib")
    file(GLOB PT_ALGLIB_FILES "${PT_DEPS_PREFIX}/lib/libalglib*.so*" "${PT_DEPS_PREFIX}/lib/libalglib*.dylib")
    install(FILES ${PT_ALGLIB_FILES} DESTINATION phasetracer COMPONENT python)
  endif()
else()
  # development build: assemble an importable package in the build tree
  set(PT_PYTHON_PACKAGE_DIR "${CMAKE_BINARY_DIR}/python/phasetracer")
  set_target_properties(_phasetracer PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${PT_PYTHON_PACKAGE_DIR}")
  # the pure-Python files are copied whenever they change, also when no C++ was rebuilt;
  # a new .py file is picked up when CMake re-runs (the list is globbed)
  file(GLOB PT_PYTHON_FILES ${CMAKE_CURRENT_SOURCE_DIR}/phasetracer/*.py)
  set(PT_PYTHON_FILES_COPIED "")
  foreach(PT_PYTHON_FILE ${PT_PYTHON_FILES})
    get_filename_component(PT_PYTHON_FILE_NAME "${PT_PYTHON_FILE}" NAME)
    add_custom_command(
      OUTPUT "${PT_PYTHON_PACKAGE_DIR}/${PT_PYTHON_FILE_NAME}"
      COMMAND ${CMAKE_COMMAND} -E copy_if_different "${PT_PYTHON_FILE}" "${PT_PYTHON_PACKAGE_DIR}/${PT_PYTHON_FILE_NAME}"
      DEPENDS "${PT_PYTHON_FILE}"
      COMMENT "Copying ${PT_PYTHON_FILE_NAME} into ${PT_PYTHON_PACKAGE_DIR}")
    list(APPEND PT_PYTHON_FILES_COPIED "${PT_PYTHON_PACKAGE_DIR}/${PT_PYTHON_FILE_NAME}")
  endforeach()
  add_custom_target(phasetracer_python_files ALL DEPENDS ${PT_PYTHON_FILES_COPIED})
  add_dependencies(_phasetracer phasetracer_python_files)
  message(STATUS "Python package will be built in ${CMAKE_BINARY_DIR}/python (add it to PYTHONPATH)")
endif()
