# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.25

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
CMAKE_COMMAND = /cvmfs/larsoft.opensciencegrid.org/products/cmake/v3_25_2/Linux64bit+3.10-2.17/bin/cmake

# The command to remove a file.
RM = /cvmfs/larsoft.opensciencegrid.org/products/cmake/v3_25_2/Linux64bit+3.10-2.17/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /exp/sbnd/app/users/theobal1/CC1muAnalysis/v09_75_03/srcs

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /exp/sbnd/app/users/theobal1/CC1muAnalysis/v09_75_03/build_slf7.x86_64

# Utility rule file for NightlyCoverage.

# Include any custom commands dependencies for this target.
include CMakeFiles/NightlyCoverage.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/NightlyCoverage.dir/progress.make

CMakeFiles/NightlyCoverage:
	/cvmfs/larsoft.opensciencegrid.org/products/cmake/v3_25_2/Linux64bit+3.10-2.17/bin/ctest -D NightlyCoverage

NightlyCoverage: CMakeFiles/NightlyCoverage
NightlyCoverage: CMakeFiles/NightlyCoverage.dir/build.make
.PHONY : NightlyCoverage

# Rule to build all files generated by this target.
CMakeFiles/NightlyCoverage.dir/build: NightlyCoverage
.PHONY : CMakeFiles/NightlyCoverage.dir/build

CMakeFiles/NightlyCoverage.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/NightlyCoverage.dir/cmake_clean.cmake
.PHONY : CMakeFiles/NightlyCoverage.dir/clean

CMakeFiles/NightlyCoverage.dir/depend:
	cd /exp/sbnd/app/users/theobal1/CC1muAnalysis/v09_75_03/build_slf7.x86_64 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /exp/sbnd/app/users/theobal1/CC1muAnalysis/v09_75_03/srcs /exp/sbnd/app/users/theobal1/CC1muAnalysis/v09_75_03/srcs /exp/sbnd/app/users/theobal1/CC1muAnalysis/v09_75_03/build_slf7.x86_64 /exp/sbnd/app/users/theobal1/CC1muAnalysis/v09_75_03/build_slf7.x86_64 /exp/sbnd/app/users/theobal1/CC1muAnalysis/v09_75_03/build_slf7.x86_64/CMakeFiles/NightlyCoverage.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/NightlyCoverage.dir/depend
