#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "meshcut::meshcut" for configuration ""
set_property(TARGET meshcut::meshcut APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(meshcut::meshcut PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libmeshcut.a"
  )

list(APPEND _cmake_import_check_targets meshcut::meshcut )
list(APPEND _cmake_import_check_files_for_meshcut::meshcut "${_IMPORT_PREFIX}/lib/libmeshcut.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
