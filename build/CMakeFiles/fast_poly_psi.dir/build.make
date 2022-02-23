# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

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
CMAKE_SOURCE_DIR = /home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/build

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/build

# Include any dependencies generated for this target.
include CMakeFiles/fast_poly_psi.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/fast_poly_psi.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/fast_poly_psi.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/fast_poly_psi.dir/flags.make

CMakeFiles/fast_poly_psi.dir/home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/fast_poly_psi.cpp.o: CMakeFiles/fast_poly_psi.dir/flags.make
CMakeFiles/fast_poly_psi.dir/home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/fast_poly_psi.cpp.o: /home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/fast_poly_psi.cpp
CMakeFiles/fast_poly_psi.dir/home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/fast_poly_psi.cpp.o: CMakeFiles/fast_poly_psi.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/fast_poly_psi.dir/home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/fast_poly_psi.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/fast_poly_psi.dir/home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/fast_poly_psi.cpp.o -MF CMakeFiles/fast_poly_psi.dir/home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/fast_poly_psi.cpp.o.d -o CMakeFiles/fast_poly_psi.dir/home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/fast_poly_psi.cpp.o -c /home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/fast_poly_psi.cpp

CMakeFiles/fast_poly_psi.dir/home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/fast_poly_psi.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fast_poly_psi.dir/home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/fast_poly_psi.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/fast_poly_psi.cpp > CMakeFiles/fast_poly_psi.dir/home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/fast_poly_psi.cpp.i

CMakeFiles/fast_poly_psi.dir/home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/fast_poly_psi.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fast_poly_psi.dir/home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/fast_poly_psi.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/fast_poly_psi.cpp -o CMakeFiles/fast_poly_psi.dir/home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/fast_poly_psi.cpp.s

# Object files for target fast_poly_psi
fast_poly_psi_OBJECTS = \
"CMakeFiles/fast_poly_psi.dir/home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/fast_poly_psi.cpp.o"

# External object files for target fast_poly_psi
fast_poly_psi_EXTERNAL_OBJECTS =

fast_poly_psi: CMakeFiles/fast_poly_psi.dir/home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/fast_poly_psi.cpp.o
fast_poly_psi: CMakeFiles/fast_poly_psi.dir/build.make
fast_poly_psi: /home/postgrad/19/mlwu/mili/SEAL/install/lib/libseal-3.6.a
fast_poly_psi: CMakeFiles/fast_poly_psi.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable fast_poly_psi"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/fast_poly_psi.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/fast_poly_psi.dir/build: fast_poly_psi
.PHONY : CMakeFiles/fast_poly_psi.dir/build

CMakeFiles/fast_poly_psi.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/fast_poly_psi.dir/cmake_clean.cmake
.PHONY : CMakeFiles/fast_poly_psi.dir/clean

CMakeFiles/fast_poly_psi.dir/depend:
	cd /home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/build /home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/build /home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/build /home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/build /home/postgrad/19/mlwu/mili/SEAL/native/examples/PolyPSI/build/CMakeFiles/fast_poly_psi.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/fast_poly_psi.dir/depend

