# Build script for the ALGLIB C++ free edition, which ships without one.
# Used by scripts/install_dependencies.sh --local: copied next to the extracted src/ directory.

cmake_minimum_required(VERSION 3.13)
project(alglib CXX)

file(GLOB ALGLIB_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
file(GLOB ALGLIB_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h)

add_library(alglib ${ALGLIB_SOURCES})
set_target_properties(alglib PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(alglib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

include(GNUInstallDirs)
install(TARGETS alglib ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${ALGLIB_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/alglib)
