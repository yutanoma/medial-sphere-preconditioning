cmake_minimum_required(VERSION 3.10.0)

project(gc_polyscope_example)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Check if Boost is installed on the system
find_package(Boost QUIET)
if(Boost_FOUND)
    message(STATUS "Boost found: ${Boost_VERSION}")
    message(STATUS "Boost include dirs: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost libraries: ${Boost_LIBRARIES}")
else()
    message(STATUS "Boost not found on system - project will use header-only alternatives")
endif()

### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Print the build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release" FORCE)
endif()
message(STATUS "cmake build type: ${CMAKE_BUILD_TYPE}")

### Configure the compiler
# This is a basic, decent setup that should do something sane on most compilers

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

  # using Clang (linux or apple) or GCC
  message("Using clang/gcc compiler flags")
  SET(BASE_CXX_FLAGS "-std=c++11 -Wall -Wextra")
  SET(DISABLED_WARNINGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations -Wno-missing-braces -Wno-unused-private-field")
  SET(TRACE_INCLUDES " -H -Wno-error=unused-command-line-argument")

  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message("Setting clang-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=3 -fcolor-diagnostics")
    SET(CMAKE_CXX_FLAGS_DEBUG          "-g3 -fsanitize=address -fno-limit-debug-info")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message("Setting gcc-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fmax-errors=5")
    SET(CMAKE_CXX_FLAGS_DEBUG          "-g3")
    SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} -Wno-maybe-uninitialized -Wno-format-zero-length -Wno-unused-but-set-parameter -Wno-unused-but-set-variable")
  endif()
  
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  SET(CMAKE_CXX_FLAGS_RELEASE        "-O3 -march=native -DNDEBUG")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
  message("Using Visual Studio compiler flags")
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /W4")
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /MP") # parallel build
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4267\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4244\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4305\"")  # ignore truncation on initialization
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")

  add_definitions(/D "_CRT_SECURE_NO_WARNINGS")
  add_definitions(-DNOMINMAX)
  add_definitions(-D_USE_MATH_DEFINES)
else()
  # unrecognized
  message( FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]" )
endif()


# == Deps
add_subdirectory(deps/geometry-central)
add_subdirectory(deps/polyscope)

# Add libIGL
# With these lines commented out, we use libIGL in "header only" mode.  Uncomment to recurse and get the full package.
# (and see the libIGL build documentation)
# message("\n\n == CMAKE recursively building libIGL\n")
# set(LIBIGL_WITH_EMBREE            OFF CACHE BOOL "Use Embree")
# set(LIBIGL_WITH_OPENGL            OFF CACHE BOOL "Use OpenGL")
# set(LIBIGL_WITH_OPENGL            OFF CACHE BOOL "Use OpenGL")
# set(LIBIGL_WITH_OPENGL_GLFW       OFF CACHE BOOL "Use GLFW"  )
# set(LIBIGL_WITH_OPENGL_GLFW_IMGUI OFF CACHE BOOL "Use ImGui" )
set(LIBIGL_RESTRICTED_TRIANGLE      ON CACHE BOOL "Use Triangle")
set(LIBIGL_COPYLEFT_CGAL            ON CACHE BOOL "Use CGAL")
# set(LIBIGL_COPYLEFT_TETGEN          ON CACHE BOOL "Use TetGen")
# note: the following line should come after setting triangle and tetgen
# set(LIBIGL_PREDICATES               ON CACHE BOOL "Use exact predicates")
add_subdirectory(deps/libigl)

# read the modules
add_subdirectory(modules)

# tinyad
add_subdirectory(deps/tinyad)

# == Build our project stuff
set(SRCS
  src/main.cpp
)

add_executable(main "${SRCS}")
target_include_directories(main PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/modules/include)
target_include_directories(main PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
# add the args.hxx project which we use for command line args
target_include_directories(main PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/polyscope/deps/args")
target_link_libraries(main PUBLIC geometry-central polyscope TinyAD modules triangle::triangle igl_copyleft::cgal)

# tinyad needs C++17 (not sure if that means this project does too)
target_compile_features(main PRIVATE cxx_std_17)

# knot untangling
add_executable(knot_untangling "src/knot_untangling.cpp")
target_include_directories(knot_untangling PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/modules/include)
target_include_directories(knot_untangling PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
target_include_directories(knot_untangling PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/polyscope/deps/args")
target_link_libraries(knot_untangling PUBLIC geometry-central polyscope TinyAD modules triangle::triangle igl_copyleft::cgal)

# tinyad needs C++17 (not sure if that means this project does too)
target_compile_features(knot_untangling PRIVATE cxx_std_17)