cmake_minimum_required(VERSION 3.14.7)

project(cpl_benchmark_bvh)

include(FetchContent)

set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_INSTALL_DOCS   OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)

FetchContent_Declare(benchmark
  URL "https://github.com/google/benchmark/archive/refs/tags/v1.7.1.zip"
  URL_HASH "SHA256=aeec52381284ec3752505a220d36096954c869da4573c2e1df3642d2f0a4aac6")

FetchContent_MakeAvailable(benchmark)

add_executable(cpl_bvh_bench
  ../../cpl_bvh.h
  cpl_bvh.c
  main.cpp)

target_link_libraries(cpl_bvh_bench PRIVATE benchmark benchmark_main)
