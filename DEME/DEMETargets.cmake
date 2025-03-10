# Generated by CMake

if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" LESS 2.8)
   message(FATAL_ERROR "CMake >= 2.8.0 required")
endif()
if(CMAKE_VERSION VERSION_LESS "2.8.3")
   message(FATAL_ERROR "CMake >= 2.8.3 required")
endif()
cmake_policy(PUSH)
cmake_policy(VERSION 2.8.3...3.24)
#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Protect against multiple inclusion, which would fail when already imported targets are added once more.
set(_cmake_targets_defined "")
set(_cmake_targets_not_defined "")
set(_cmake_expected_targets "")
foreach(_cmake_expected_target IN ITEMS DEME::simulator_multi_gpu DEME::DEMERuntimeDataHelper DEME::DEMERuntimeDataHelper_install)
  list(APPEND _cmake_expected_targets "${_cmake_expected_target}")
  if(TARGET "${_cmake_expected_target}")
    list(APPEND _cmake_targets_defined "${_cmake_expected_target}")
  else()
    list(APPEND _cmake_targets_not_defined "${_cmake_expected_target}")
  endif()
endforeach()
unset(_cmake_expected_target)
if(_cmake_targets_defined STREQUAL _cmake_expected_targets)
  unset(_cmake_targets_defined)
  unset(_cmake_targets_not_defined)
  unset(_cmake_expected_targets)
  unset(CMAKE_IMPORT_FILE_VERSION)
  cmake_policy(POP)
  return()
endif()
if(NOT _cmake_targets_defined STREQUAL "")
  string(REPLACE ";" ", " _cmake_targets_defined_text "${_cmake_targets_defined}")
  string(REPLACE ";" ", " _cmake_targets_not_defined_text "${_cmake_targets_not_defined}")
  message(FATAL_ERROR "Some (but not all) targets in this export set were already defined.\nTargets Defined: ${_cmake_targets_defined_text}\nTargets not yet defined: ${_cmake_targets_not_defined_text}\n")
endif()
unset(_cmake_targets_defined)
unset(_cmake_targets_not_defined)
unset(_cmake_expected_targets)


# Create imported target DEME::simulator_multi_gpu
add_library(DEME::simulator_multi_gpu STATIC IMPORTED)

set_target_properties(DEME::simulator_multi_gpu PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "/home/lukasg/dem/DEM-Engine/src;/home/lukasg/dem/build/src;/home/lukasg/dem/DEM-Engine/thirdparty/nvidia_helper_math"
  INTERFACE_LINK_LIBRARIES "CUDA::cudart;CUDA::nvrtc;CUDA::cuda_driver;DEME::DEMERuntimeDataHelper"
)

# Create imported target DEME::DEMERuntimeDataHelper
add_library(DEME::DEMERuntimeDataHelper SHARED IMPORTED)

set_target_properties(DEME::DEMERuntimeDataHelper PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "/home/lukasg/dem/DEM-Engine/src"
  INTERFACE_SOURCES "/home/lukasg/dem/DEM-Engine/src/core/utils/RuntimeData.h"
)

# Create imported target DEME::DEMERuntimeDataHelper_install
add_library(DEME::DEMERuntimeDataHelper_install SHARED IMPORTED)

set_target_properties(DEME::DEMERuntimeDataHelper_install PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "/home/lukasg/dem/DEM-Engine/src"
  INTERFACE_LINK_OPTIONS "LINKER:-soname,libDEMERuntimeDataHelper.so"
  INTERFACE_SOURCES "/home/lukasg/dem/DEM-Engine/src/core/utils/RuntimeData.h"
)

# Import target "DEME::simulator_multi_gpu" for configuration "Release"
set_property(TARGET DEME::simulator_multi_gpu APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(DEME::simulator_multi_gpu PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CUDA;CXX"
  IMPORTED_LOCATION_RELEASE "/home/lukasg/dem/build/libsimulator_multi_gpu.a"
  )

# Import target "DEME::DEMERuntimeDataHelper" for configuration "Release"
set_property(TARGET DEME::DEMERuntimeDataHelper APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(DEME::DEMERuntimeDataHelper PROPERTIES
  IMPORTED_LOCATION_RELEASE "/home/lukasg/dem/build/src/core/libDEMERuntimeDataHelper.so"
  IMPORTED_SONAME_RELEASE "libDEMERuntimeDataHelper.so"
  )

# Import target "DEME::DEMERuntimeDataHelper_install" for configuration "Release"
set_property(TARGET DEME::DEMERuntimeDataHelper_install APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(DEME::DEMERuntimeDataHelper_install PROPERTIES
  IMPORTED_LOCATION_RELEASE "/home/lukasg/dem/build/lib_install/libDEMERuntimeDataHelper.so"
  IMPORTED_SONAME_RELEASE "libDEMERuntimeDataHelper.so"
  )

# This file does not depend on other imported targets which have
# been exported from the same project but in a separate export set.

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
cmake_policy(POP)
