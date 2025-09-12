cmake_minimum_required(VERSION 3.12)
project(meshcut VERSION 1.0.0)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add library
add_library(meshcut STATIC
    src/meshcut.cpp
    src/clipper.engine.cpp
    src/clipper.offset.cpp
    src/clipper.rectclip.cpp
)

# Include directories
target_include_directories(meshcut PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Compiler-specific options
target_compile_features(meshcut PUBLIC cxx_std_17)

# Install
include(GNUInstallDirs)

install(TARGETS meshcut
    EXPORT meshcutTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT meshcutTargets
    FILE meshcutTargets.cmake
    NAMESPACE meshcut::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/meshcut
)

# Create config file
include(CMakePackageConfigHelpers)
configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/meshcutConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/meshcutConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/meshcut
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/meshcutConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/meshcutConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/meshcutConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/meshcut
)

# Examples (optional)
option(MESHCUT_BUILD_EXAMPLES "Build MeshCut examples" OFF)
if(MESHCUT_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()