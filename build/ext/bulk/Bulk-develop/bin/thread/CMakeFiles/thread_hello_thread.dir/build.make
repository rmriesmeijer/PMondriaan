# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

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
CMAKE_COMMAND = /usr/bin/cmake.exe

# The command to remove a file.
RM = /usr/bin/cmake.exe -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /cygdrive/c/users/sarit/Documents/MasterThesis/programming

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /cygdrive/c/users/sarit/Documents/MasterThesis/programming/build

# Include any dependencies generated for this target.
include ext/bulk/Bulk-develop/bin/thread/CMakeFiles/thread_hello_thread.dir/depend.make

# Include the progress variables for this target.
include ext/bulk/Bulk-develop/bin/thread/CMakeFiles/thread_hello_thread.dir/progress.make

# Include the compile flags for this target's objects.
include ext/bulk/Bulk-develop/bin/thread/CMakeFiles/thread_hello_thread.dir/flags.make

ext/bulk/Bulk-develop/bin/thread/CMakeFiles/thread_hello_thread.dir/__/examples/hello_thread.cpp.o: ext/bulk/Bulk-develop/bin/thread/CMakeFiles/thread_hello_thread.dir/flags.make
ext/bulk/Bulk-develop/bin/thread/CMakeFiles/thread_hello_thread.dir/__/examples/hello_thread.cpp.o: ../ext/bulk/Bulk-develop/backends/thread/examples/hello_thread.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cygdrive/c/users/sarit/Documents/MasterThesis/programming/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object ext/bulk/Bulk-develop/bin/thread/CMakeFiles/thread_hello_thread.dir/__/examples/hello_thread.cpp.o"
	cd /cygdrive/c/users/sarit/Documents/MasterThesis/programming/build/ext/bulk/Bulk-develop/bin/thread && /usr/bin/c++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/thread_hello_thread.dir/__/examples/hello_thread.cpp.o -c /cygdrive/c/users/sarit/Documents/MasterThesis/programming/ext/bulk/Bulk-develop/backends/thread/examples/hello_thread.cpp

ext/bulk/Bulk-develop/bin/thread/CMakeFiles/thread_hello_thread.dir/__/examples/hello_thread.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/thread_hello_thread.dir/__/examples/hello_thread.cpp.i"
	cd /cygdrive/c/users/sarit/Documents/MasterThesis/programming/build/ext/bulk/Bulk-develop/bin/thread && /usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cygdrive/c/users/sarit/Documents/MasterThesis/programming/ext/bulk/Bulk-develop/backends/thread/examples/hello_thread.cpp > CMakeFiles/thread_hello_thread.dir/__/examples/hello_thread.cpp.i

ext/bulk/Bulk-develop/bin/thread/CMakeFiles/thread_hello_thread.dir/__/examples/hello_thread.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/thread_hello_thread.dir/__/examples/hello_thread.cpp.s"
	cd /cygdrive/c/users/sarit/Documents/MasterThesis/programming/build/ext/bulk/Bulk-develop/bin/thread && /usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cygdrive/c/users/sarit/Documents/MasterThesis/programming/ext/bulk/Bulk-develop/backends/thread/examples/hello_thread.cpp -o CMakeFiles/thread_hello_thread.dir/__/examples/hello_thread.cpp.s

# Object files for target thread_hello_thread
thread_hello_thread_OBJECTS = \
"CMakeFiles/thread_hello_thread.dir/__/examples/hello_thread.cpp.o"

# External object files for target thread_hello_thread
thread_hello_thread_EXTERNAL_OBJECTS =

ext/bulk/Bulk-develop/bin/thread/thread_hello_thread.exe: ext/bulk/Bulk-develop/bin/thread/CMakeFiles/thread_hello_thread.dir/__/examples/hello_thread.cpp.o
ext/bulk/Bulk-develop/bin/thread/thread_hello_thread.exe: ext/bulk/Bulk-develop/bin/thread/CMakeFiles/thread_hello_thread.dir/build.make
ext/bulk/Bulk-develop/bin/thread/thread_hello_thread.exe: ext/bulk/Bulk-develop/bin/thread/CMakeFiles/thread_hello_thread.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/cygdrive/c/users/sarit/Documents/MasterThesis/programming/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable thread_hello_thread.exe"
	cd /cygdrive/c/users/sarit/Documents/MasterThesis/programming/build/ext/bulk/Bulk-develop/bin/thread && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/thread_hello_thread.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
ext/bulk/Bulk-develop/bin/thread/CMakeFiles/thread_hello_thread.dir/build: ext/bulk/Bulk-develop/bin/thread/thread_hello_thread.exe

.PHONY : ext/bulk/Bulk-develop/bin/thread/CMakeFiles/thread_hello_thread.dir/build

ext/bulk/Bulk-develop/bin/thread/CMakeFiles/thread_hello_thread.dir/clean:
	cd /cygdrive/c/users/sarit/Documents/MasterThesis/programming/build/ext/bulk/Bulk-develop/bin/thread && $(CMAKE_COMMAND) -P CMakeFiles/thread_hello_thread.dir/cmake_clean.cmake
.PHONY : ext/bulk/Bulk-develop/bin/thread/CMakeFiles/thread_hello_thread.dir/clean

ext/bulk/Bulk-develop/bin/thread/CMakeFiles/thread_hello_thread.dir/depend:
	cd /cygdrive/c/users/sarit/Documents/MasterThesis/programming/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /cygdrive/c/users/sarit/Documents/MasterThesis/programming /cygdrive/c/users/sarit/Documents/MasterThesis/programming/ext/bulk/Bulk-develop/backends/thread/build /cygdrive/c/users/sarit/Documents/MasterThesis/programming/build /cygdrive/c/users/sarit/Documents/MasterThesis/programming/build/ext/bulk/Bulk-develop/bin/thread /cygdrive/c/users/sarit/Documents/MasterThesis/programming/build/ext/bulk/Bulk-develop/bin/thread/CMakeFiles/thread_hello_thread.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : ext/bulk/Bulk-develop/bin/thread/CMakeFiles/thread_hello_thread.dir/depend

