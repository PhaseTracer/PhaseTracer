# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

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
CMAKE_COMMAND = /home/Alfie/.local/lib/python2.7/site-packages/cmake/data/bin/cmake

# The command to remove a file.
RM = /home/Alfie/.local/lib/python2.7/site-packages/cmake/data/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/e/Software/PhaseTracer

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/e/Software/PhaseTracer/build

# Include any dependencies generated for this target.
include src/CMakeFiles/phasetracer.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/phasetracer.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/phasetracer.dir/flags.make

src/CMakeFiles/phasetracer.dir/h_bar_expansion.cpp.o: src/CMakeFiles/phasetracer.dir/flags.make
src/CMakeFiles/phasetracer.dir/h_bar_expansion.cpp.o: ../src/h_bar_expansion.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/e/Software/PhaseTracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/phasetracer.dir/h_bar_expansion.cpp.o"
	cd /mnt/e/Software/PhaseTracer/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/phasetracer.dir/h_bar_expansion.cpp.o -c /mnt/e/Software/PhaseTracer/src/h_bar_expansion.cpp

src/CMakeFiles/phasetracer.dir/h_bar_expansion.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/phasetracer.dir/h_bar_expansion.cpp.i"
	cd /mnt/e/Software/PhaseTracer/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/e/Software/PhaseTracer/src/h_bar_expansion.cpp > CMakeFiles/phasetracer.dir/h_bar_expansion.cpp.i

src/CMakeFiles/phasetracer.dir/h_bar_expansion.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/phasetracer.dir/h_bar_expansion.cpp.s"
	cd /mnt/e/Software/PhaseTracer/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/e/Software/PhaseTracer/src/h_bar_expansion.cpp -o CMakeFiles/phasetracer.dir/h_bar_expansion.cpp.s

src/CMakeFiles/phasetracer.dir/phase_finder.cpp.o: src/CMakeFiles/phasetracer.dir/flags.make
src/CMakeFiles/phasetracer.dir/phase_finder.cpp.o: ../src/phase_finder.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/e/Software/PhaseTracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/CMakeFiles/phasetracer.dir/phase_finder.cpp.o"
	cd /mnt/e/Software/PhaseTracer/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/phasetracer.dir/phase_finder.cpp.o -c /mnt/e/Software/PhaseTracer/src/phase_finder.cpp

src/CMakeFiles/phasetracer.dir/phase_finder.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/phasetracer.dir/phase_finder.cpp.i"
	cd /mnt/e/Software/PhaseTracer/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/e/Software/PhaseTracer/src/phase_finder.cpp > CMakeFiles/phasetracer.dir/phase_finder.cpp.i

src/CMakeFiles/phasetracer.dir/phase_finder.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/phasetracer.dir/phase_finder.cpp.s"
	cd /mnt/e/Software/PhaseTracer/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/e/Software/PhaseTracer/src/phase_finder.cpp -o CMakeFiles/phasetracer.dir/phase_finder.cpp.s

src/CMakeFiles/phasetracer.dir/transition_finder.cpp.o: src/CMakeFiles/phasetracer.dir/flags.make
src/CMakeFiles/phasetracer.dir/transition_finder.cpp.o: ../src/transition_finder.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/e/Software/PhaseTracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/CMakeFiles/phasetracer.dir/transition_finder.cpp.o"
	cd /mnt/e/Software/PhaseTracer/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/phasetracer.dir/transition_finder.cpp.o -c /mnt/e/Software/PhaseTracer/src/transition_finder.cpp

src/CMakeFiles/phasetracer.dir/transition_finder.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/phasetracer.dir/transition_finder.cpp.i"
	cd /mnt/e/Software/PhaseTracer/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/e/Software/PhaseTracer/src/transition_finder.cpp > CMakeFiles/phasetracer.dir/transition_finder.cpp.i

src/CMakeFiles/phasetracer.dir/transition_finder.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/phasetracer.dir/transition_finder.cpp.s"
	cd /mnt/e/Software/PhaseTracer/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/e/Software/PhaseTracer/src/transition_finder.cpp -o CMakeFiles/phasetracer.dir/transition_finder.cpp.s

# Object files for target phasetracer
phasetracer_OBJECTS = \
"CMakeFiles/phasetracer.dir/h_bar_expansion.cpp.o" \
"CMakeFiles/phasetracer.dir/phase_finder.cpp.o" \
"CMakeFiles/phasetracer.dir/transition_finder.cpp.o"

# External object files for target phasetracer
phasetracer_EXTERNAL_OBJECTS =

../lib/libphasetracer.so: src/CMakeFiles/phasetracer.dir/h_bar_expansion.cpp.o
../lib/libphasetracer.so: src/CMakeFiles/phasetracer.dir/phase_finder.cpp.o
../lib/libphasetracer.so: src/CMakeFiles/phasetracer.dir/transition_finder.cpp.o
../lib/libphasetracer.so: src/CMakeFiles/phasetracer.dir/build.make
../lib/libphasetracer.so: /usr/lib/x86_64-linux-gnu/libnlopt.so
../lib/libphasetracer.so: /usr/lib/x86_64-linux-gnu/libboost_filesystem.so
../lib/libphasetracer.so: /usr/lib/x86_64-linux-gnu/libboost_log.so
../lib/libphasetracer.so: /usr/lib/x86_64-linux-gnu/libboost_system.so
../lib/libphasetracer.so: /usr/lib/x86_64-linux-gnu/libboost_log_setup.so
../lib/libphasetracer.so: /usr/lib/x86_64-linux-gnu/libboost_date_time.so
../lib/libphasetracer.so: /usr/lib/x86_64-linux-gnu/libboost_thread.so
../lib/libphasetracer.so: /usr/lib/x86_64-linux-gnu/libboost_regex.so
../lib/libphasetracer.so: /usr/lib/x86_64-linux-gnu/libboost_chrono.so
../lib/libphasetracer.so: /usr/lib/x86_64-linux-gnu/libboost_atomic.so
../lib/libphasetracer.so: ../EffectivePotential/lib/libeffectivepotential.so
../lib/libphasetracer.so: /usr/lib/x86_64-linux-gnu/libboost_regex.so
../lib/libphasetracer.so: /usr/lib/x86_64-linux-gnu/libboost_chrono.so
../lib/libphasetracer.so: /usr/lib/x86_64-linux-gnu/libboost_atomic.so
../lib/libphasetracer.so: /usr/lib/x86_64-linux-gnu/libalglib.so
../lib/libphasetracer.so: src/CMakeFiles/phasetracer.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/e/Software/PhaseTracer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX shared library ../../lib/libphasetracer.so"
	cd /mnt/e/Software/PhaseTracer/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/phasetracer.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/phasetracer.dir/build: ../lib/libphasetracer.so

.PHONY : src/CMakeFiles/phasetracer.dir/build

src/CMakeFiles/phasetracer.dir/clean:
	cd /mnt/e/Software/PhaseTracer/build/src && $(CMAKE_COMMAND) -P CMakeFiles/phasetracer.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/phasetracer.dir/clean

src/CMakeFiles/phasetracer.dir/depend:
	cd /mnt/e/Software/PhaseTracer/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/e/Software/PhaseTracer /mnt/e/Software/PhaseTracer/src /mnt/e/Software/PhaseTracer/build /mnt/e/Software/PhaseTracer/build/src /mnt/e/Software/PhaseTracer/build/src/CMakeFiles/phasetracer.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/phasetracer.dir/depend

