# Copyright (c) 2021, SBEL GPU Development Team
# Copyright (c) 2021, University of Wisconsin - Madison
# 
#	SPDX-License-Identifier: BSD-3-Clause


# Installable DEMEConfig.cmake (base file) 
#
# Provides the following package variables:
#
# DEME_INCLUDE_DIRS
# DEME_LIBRARIES
#
# Requires the following packages to be loaded first:
#
# CUDAToolkit
#


if (NOT CUDAToolkit_FOUND AND NOT DEME_IGNORE_MISSING_CUDA)
	message(SEND_ERROR "DEME requires the CUDAToolkit package in order to build. To suppress this error, set the variable DEME_IGNORE_MISSING_CUDA.")
endif()

cmake_path(GET CMAKE_CURRENT_LIST_FILE PARENT_PATH DEMECMakeDir)

if (NOT TARGET simulator_multi_gpu AND NOT DEME_BINARY_DIR)
	include("${DEMECMakeDir}/DEMETargets.cmake")
endif()

set(DEME_INCLUDE_DIRS ${DEMECMakeDir}/../../../include)
set(DEME_DATA_DIRS ${DEMECMakeDir}/../../../include/../share/DEME/data)
set(DEME_LIBRARIES DEME::simulator_multi_gpu)

# The info on whether it is compiled with ChPF on
set(DEME_WITH_CHPF "OFF")

