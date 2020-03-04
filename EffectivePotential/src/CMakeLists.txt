# Create the library target

file(GLOB SRC_FILES *.cpp)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
link_libraries(${ALGLIB_LIBRARIES})
add_library(effectivepotential SHARED ${SRC_FILES})

# add header files here so they are available in the cmake config

target_include_directories(effectivepotential INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/effectivepotential>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
)

# helpers for making config

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# package information

write_basic_package_version_file(
    effectivepotentialConfigVersion.cmake
    VERSION ${PACKAGE_VERSION}
    COMPATIBILITY AnyNewerVersion
    )

# export config file in package registry (usually ~/.cmake/packages)

install(TARGETS effectivepotential
        EXPORT effectivepotential-config
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(TARGETS effectivepotential FILE effectivepotential-config.cmake)
export(PACKAGE EffectivePotential)
