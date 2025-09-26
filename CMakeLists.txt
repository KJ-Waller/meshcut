cmake_minimum_required(VERSION 3.16)
project(meshcut VERSION 1.0.0 LANGUAGES CXX)

# Set C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(include)
include_directories(third_party)

# Add mapbox namespace include for earcut
add_definitions(-DMAPBOX_EARCUT_HEADER_ONLY)

# Compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -Wall")

# Default to Release build
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# MeshCut library (header-only)
add_library(meshcut INTERFACE)

target_include_directories(meshcut INTERFACE include)
target_include_directories(meshcut INTERFACE third_party)

# Test executables
add_executable(meshcut_test
    test/test_basic.cpp
)

add_executable(meshcut_advanced_test
    test/test_advanced.cpp
)

add_executable(meshcut_visualization_test
    test/test_visualization.cpp
)

add_executable(meshcut_performance_test
    test/test_performance.cpp
)

add_executable(meshcut_spatial_comparison_test
    test/test_spatial_comparison.cpp
)

add_executable(meshcut_visualization_complex_test
    test/test_visualization_complex.cpp
)

# Find curl and zlib for real vector tile tests
find_package(PkgConfig REQUIRED)
pkg_check_modules(CURL REQUIRED libcurl)
find_package(ZLIB REQUIRED)

add_executable(meshcut_real_polygon_extraction_test
    test/test_real_polygon_extraction.cpp
)

target_link_libraries(meshcut_test meshcut)
target_link_libraries(meshcut_advanced_test meshcut)
target_link_libraries(meshcut_visualization_test meshcut)
target_link_libraries(meshcut_performance_test meshcut)
target_link_libraries(meshcut_spatial_comparison_test meshcut)
target_link_libraries(meshcut_visualization_complex_test meshcut)
target_link_libraries(meshcut_real_polygon_extraction_test meshcut ${CURL_LIBRARIES} ${ZLIB_LIBRARIES})

# Include curl headers for real polygon extraction test
target_include_directories(meshcut_real_polygon_extraction_test PRIVATE ${CURL_INCLUDE_DIRS})
target_compile_options(meshcut_real_polygon_extraction_test PRIVATE ${CURL_CFLAGS_OTHER})

# Enable testing
enable_testing()
add_test(NAME basic_test COMMAND meshcut_test)
add_test(NAME advanced_test COMMAND meshcut_advanced_test)
add_test(NAME visualization_test COMMAND meshcut_visualization_test)
add_test(NAME performance_test COMMAND meshcut_performance_test)