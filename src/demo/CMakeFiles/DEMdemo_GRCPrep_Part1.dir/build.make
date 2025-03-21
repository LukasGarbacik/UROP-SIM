# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/lukasg/dem/DEM-Engine

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/lukasg/dem/build

# Include any dependencies generated for this target.
include src/demo/CMakeFiles/DEMdemo_GRCPrep_Part1.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/demo/CMakeFiles/DEMdemo_GRCPrep_Part1.dir/compiler_depend.make

# Include the progress variables for this target.
include src/demo/CMakeFiles/DEMdemo_GRCPrep_Part1.dir/progress.make

# Include the compile flags for this target's objects.
include src/demo/CMakeFiles/DEMdemo_GRCPrep_Part1.dir/flags.make

src/demo/CMakeFiles/DEMdemo_GRCPrep_Part1.dir/DEMdemo_GRCPrep_Part1.cpp.o: src/demo/CMakeFiles/DEMdemo_GRCPrep_Part1.dir/flags.make
src/demo/CMakeFiles/DEMdemo_GRCPrep_Part1.dir/DEMdemo_GRCPrep_Part1.cpp.o: /home/lukasg/dem/DEM-Engine/src/demo/DEMdemo_GRCPrep_Part1.cpp
src/demo/CMakeFiles/DEMdemo_GRCPrep_Part1.dir/DEMdemo_GRCPrep_Part1.cpp.o: src/demo/CMakeFiles/DEMdemo_GRCPrep_Part1.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lukasg/dem/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/demo/CMakeFiles/DEMdemo_GRCPrep_Part1.dir/DEMdemo_GRCPrep_Part1.cpp.o"
	cd /home/lukasg/dem/build/src/demo && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/demo/CMakeFiles/DEMdemo_GRCPrep_Part1.dir/DEMdemo_GRCPrep_Part1.cpp.o -MF CMakeFiles/DEMdemo_GRCPrep_Part1.dir/DEMdemo_GRCPrep_Part1.cpp.o.d -o CMakeFiles/DEMdemo_GRCPrep_Part1.dir/DEMdemo_GRCPrep_Part1.cpp.o -c /home/lukasg/dem/DEM-Engine/src/demo/DEMdemo_GRCPrep_Part1.cpp

src/demo/CMakeFiles/DEMdemo_GRCPrep_Part1.dir/DEMdemo_GRCPrep_Part1.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DEMdemo_GRCPrep_Part1.dir/DEMdemo_GRCPrep_Part1.cpp.i"
	cd /home/lukasg/dem/build/src/demo && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lukasg/dem/DEM-Engine/src/demo/DEMdemo_GRCPrep_Part1.cpp > CMakeFiles/DEMdemo_GRCPrep_Part1.dir/DEMdemo_GRCPrep_Part1.cpp.i

src/demo/CMakeFiles/DEMdemo_GRCPrep_Part1.dir/DEMdemo_GRCPrep_Part1.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DEMdemo_GRCPrep_Part1.dir/DEMdemo_GRCPrep_Part1.cpp.s"
	cd /home/lukasg/dem/build/src/demo && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lukasg/dem/DEM-Engine/src/demo/DEMdemo_GRCPrep_Part1.cpp -o CMakeFiles/DEMdemo_GRCPrep_Part1.dir/DEMdemo_GRCPrep_Part1.cpp.s

# Object files for target DEMdemo_GRCPrep_Part1
DEMdemo_GRCPrep_Part1_OBJECTS = \
"CMakeFiles/DEMdemo_GRCPrep_Part1.dir/DEMdemo_GRCPrep_Part1.cpp.o"

# External object files for target DEMdemo_GRCPrep_Part1
DEMdemo_GRCPrep_Part1_EXTERNAL_OBJECTS =

bin/DEMdemo_GRCPrep_Part1: src/demo/CMakeFiles/DEMdemo_GRCPrep_Part1.dir/DEMdemo_GRCPrep_Part1.cpp.o
bin/DEMdemo_GRCPrep_Part1: src/demo/CMakeFiles/DEMdemo_GRCPrep_Part1.dir/build.make
bin/DEMdemo_GRCPrep_Part1: libsimulator_multi_gpu.a
bin/DEMdemo_GRCPrep_Part1: /usr/local/cuda-12.1/targets/x86_64-linux/lib/libcudart.so
bin/DEMdemo_GRCPrep_Part1: /usr/local/cuda-12.1/targets/x86_64-linux/lib/libnvrtc.so
bin/DEMdemo_GRCPrep_Part1: /usr/local/cuda-12.1/targets/x86_64-linux/lib/libnvJitLink.so
bin/DEMdemo_GRCPrep_Part1: /usr/lib64/libcuda.so
bin/DEMdemo_GRCPrep_Part1: src/core/libDEMERuntimeDataHelper.so
bin/DEMdemo_GRCPrep_Part1: src/demo/CMakeFiles/DEMdemo_GRCPrep_Part1.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/lukasg/dem/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../bin/DEMdemo_GRCPrep_Part1"
	cd /home/lukasg/dem/build/src/demo && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/DEMdemo_GRCPrep_Part1.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/demo/CMakeFiles/DEMdemo_GRCPrep_Part1.dir/build: bin/DEMdemo_GRCPrep_Part1
.PHONY : src/demo/CMakeFiles/DEMdemo_GRCPrep_Part1.dir/build

src/demo/CMakeFiles/DEMdemo_GRCPrep_Part1.dir/clean:
	cd /home/lukasg/dem/build/src/demo && $(CMAKE_COMMAND) -P CMakeFiles/DEMdemo_GRCPrep_Part1.dir/cmake_clean.cmake
.PHONY : src/demo/CMakeFiles/DEMdemo_GRCPrep_Part1.dir/clean

src/demo/CMakeFiles/DEMdemo_GRCPrep_Part1.dir/depend:
	cd /home/lukasg/dem/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/lukasg/dem/DEM-Engine /home/lukasg/dem/DEM-Engine/src/demo /home/lukasg/dem/build /home/lukasg/dem/build/src/demo /home/lukasg/dem/build/src/demo/CMakeFiles/DEMdemo_GRCPrep_Part1.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/demo/CMakeFiles/DEMdemo_GRCPrep_Part1.dir/depend

