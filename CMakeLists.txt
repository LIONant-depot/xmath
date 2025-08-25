cmake_minimum_required(VERSION 3.10)
project(xmath_unit_test LANGUAGES CXX)
set(TARGET_PROJECT "xmath_unit_test")
set(COMPONENT_REGISTRY "")


# disable strict string conformance. AMD compresator requires it!
if(MSVC)
  add_compile_options(/Zc:strictStrings-)
endif()

include(FetchContent)
FetchContent_Declare(
  xcmake_tools
  GIT_REPOSITORY https://github.com/LIONant-depot/xcmake_tools.git
  GIT_TAG main
  GIT_SHALLOW TRUE
  SOURCE_DIR "${CMAKE_BINARY_DIR}/_deps/xcmake_tools"
)
FetchContent_MakeAvailable(xcmake_tools)
include(${CMAKE_BINARY_DIR}/_deps/xcmake_tools/Common.cmake)

add_executable(${TARGET_PROJECT}
  "source/unittests/xmath_fmat3_unittest.h"
  "source/unittests/xmath_fmat4_unittest.h"
  "source/unittests/xmath_fquat_unittest.h"
  "source/unittests/xmath_fvec2_unittest.h"
  "source/unittests/xmath_fvec3_unittest.h"
  "source/unittests/xmath_fvec4_unittest.h"
  "source/unittests/xmath_radian3_unittest.h"
  "source/unittests/xmath_unittest.h"
  "source/unittests/main.cpp"
)

source_group("Unit Tests" FILES
  "source/unittests/xmath_fmat3_unittest.h"
  "source/unittests/xmath_fmat4_unittest.h"
  "source/unittests/xmath_fquat_unittest.h"
  "source/unittests/xmath_fvec2_unittest.h"
  "source/unittests/xmath_fvec3_unittest.h"
  "source/unittests/xmath_fvec4_unittest.h"
  "source/unittests/xmath_radian3_unittest.h"
  "source/unittests/xmath_unittest.h"
)

source_group("" FILES
  "source/unittests/main.cpp"
)

add_subdirectory("build/dependency" "${CMAKE_CURRENT_BINARY_DIR}/xmath_unit_test")

ProcessComponents()