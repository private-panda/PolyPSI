# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
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
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/build

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/build

# Include any dependencies generated for this target.
include CMakeFiles/low_poly_psi.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/low_poly_psi.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/low_poly_psi.dir/flags.make

CMakeFiles/low_poly_psi.dir/home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/low_poly_psi.cpp.o: CMakeFiles/low_poly_psi.dir/flags.make
CMakeFiles/low_poly_psi.dir/home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/low_poly_psi.cpp.o: /home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/low_poly_psi.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/low_poly_psi.dir/home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/low_poly_psi.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/low_poly_psi.dir/home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/low_poly_psi.cpp.o -c /home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/low_poly_psi.cpp

CMakeFiles/low_poly_psi.dir/home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/low_poly_psi.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/low_poly_psi.dir/home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/low_poly_psi.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/low_poly_psi.cpp > CMakeFiles/low_poly_psi.dir/home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/low_poly_psi.cpp.i

CMakeFiles/low_poly_psi.dir/home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/low_poly_psi.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/low_poly_psi.dir/home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/low_poly_psi.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/low_poly_psi.cpp -o CMakeFiles/low_poly_psi.dir/home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/low_poly_psi.cpp.s

# Object files for target low_poly_psi
low_poly_psi_OBJECTS = \
"CMakeFiles/low_poly_psi.dir/home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/low_poly_psi.cpp.o"

# External object files for target low_poly_psi
low_poly_psi_EXTERNAL_OBJECTS =

low_poly_psi: CMakeFiles/low_poly_psi.dir/home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/low_poly_psi.cpp.o
low_poly_psi: CMakeFiles/low_poly_psi.dir/build.make
low_poly_psi: /usr/local/lib/libseal-3.6.a
low_poly_psi: CMakeFiles/low_poly_psi.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable low_poly_psi"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/low_poly_psi.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/low_poly_psi.dir/build: low_poly_psi

.PHONY : CMakeFiles/low_poly_psi.dir/build

CMakeFiles/low_poly_psi.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/low_poly_psi.dir/cmake_clean.cmake
.PHONY : CMakeFiles/low_poly_psi.dir/clean

CMakeFiles/low_poly_psi.dir/depend:
	cd /home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/build /home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/build /home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/build /home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/build /home/mingli/Documents/SEAL/native/examples/PolyPSI-OPRF/build/CMakeFiles/low_poly_psi.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/low_poly_psi.dir/depend

