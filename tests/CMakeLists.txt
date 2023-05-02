cmake_minimum_required(VERSION 3.14.7)

project(cpl_tests)

include(FetchContent)

FetchContent_Declare(googletest
  URL "https://github.com/google/googletest/archive/refs/tags/v1.13.0.zip"
  URL_HASH "SHA256=ffa17fbc5953900994e2deec164bb8949879ea09b411e07f215bfbb1f87f4632")

FetchContent_MakeAvailable(googletest)

add_library(cpl_testlib
  ../cpl_math.h
  cpl_math.c)

if(CMAKE_COMPILER_IS_GNUCC)
  target_compile_options(cpl_testlib
    PRIVATE
      -Wall -Wextra -Werror -Wfatal-errors -Wconversion -Wstrict-aliasing -std=c89 -pedantic -ansi -fanalyzer)
endif(CMAKE_COMPILER_IS_GNUCC)

add_executable(cpl_tests
  math_tests.cpp)

target_link_libraries(cpl_tests PRIVATE gtest gtest_main cpl_testlib)

if(CMAKE_COMPILER_IS_GNUCXX)
  target_compile_options(cpl_tests PRIVATE -Wall -Wextra -Werror -Wfatal-errors)
endif(CMAKE_COMPILER_IS_GNUCXX)
