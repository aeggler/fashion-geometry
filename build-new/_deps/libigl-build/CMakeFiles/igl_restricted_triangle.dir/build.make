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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.26.3/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.26.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/annaeggler/Desktop/Documents/Studium/Masterarbeit/fashion-descriptors

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/annaeggler/Desktop/Documents/Studium/Masterarbeit/fashion-descriptors/build-new

# Utility rule file for igl_restricted_triangle.

# Include any custom commands dependencies for this target.
include _deps/libigl-build/CMakeFiles/igl_restricted_triangle.dir/compiler_depend.make

# Include the progress variables for this target.
include _deps/libigl-build/CMakeFiles/igl_restricted_triangle.dir/progress.make

igl_restricted_triangle: _deps/libigl-build/CMakeFiles/igl_restricted_triangle.dir/build.make
.PHONY : igl_restricted_triangle

# Rule to build all files generated by this target.
_deps/libigl-build/CMakeFiles/igl_restricted_triangle.dir/build: igl_restricted_triangle
.PHONY : _deps/libigl-build/CMakeFiles/igl_restricted_triangle.dir/build

_deps/libigl-build/CMakeFiles/igl_restricted_triangle.dir/clean:
	cd /Users/annaeggler/Desktop/Documents/Studium/Masterarbeit/fashion-descriptors/build-new/_deps/libigl-build && $(CMAKE_COMMAND) -P CMakeFiles/igl_restricted_triangle.dir/cmake_clean.cmake
.PHONY : _deps/libigl-build/CMakeFiles/igl_restricted_triangle.dir/clean

_deps/libigl-build/CMakeFiles/igl_restricted_triangle.dir/depend:
	cd /Users/annaeggler/Desktop/Documents/Studium/Masterarbeit/fashion-descriptors/build-new && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/annaeggler/Desktop/Documents/Studium/Masterarbeit/fashion-descriptors /Users/annaeggler/Desktop/Documents/Studium/Masterarbeit/fashion-descriptors/build-new/_deps/libigl-src /Users/annaeggler/Desktop/Documents/Studium/Masterarbeit/fashion-descriptors/build-new /Users/annaeggler/Desktop/Documents/Studium/Masterarbeit/fashion-descriptors/build-new/_deps/libigl-build /Users/annaeggler/Desktop/Documents/Studium/Masterarbeit/fashion-descriptors/build-new/_deps/libigl-build/CMakeFiles/igl_restricted_triangle.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : _deps/libigl-build/CMakeFiles/igl_restricted_triangle.dir/depend

