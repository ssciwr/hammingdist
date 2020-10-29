project(covid-tda LANGUAGES CXX)

cmake_minimum_required(VERSION 3.11)

set(HAMMING_WITH_OPENMP
    no
    CACHE BOOL "Build with OpenMP support (multithreading)")

# Build library
add_library(hamming STATIC hamming.cc)
target_compile_options(hamming PUBLIC -fPIC)
if(HAMMING_WITH_OPENMP)
  find_package(OpenMP REQUIRED)
  target_compile_definitions(hamming PUBLIC HAMMING_WITH_OPENMP)
  target_link_libraries(hamming PRIVATE OpenMP::OpenMP_CXX)
endif()

# Build library benchmarks
set(BENCHMARK_ENABLE_GTEST_TESTS OFF)
set(BENCHMARK_ENABLE_TESTING OFF)
add_subdirectory(benchmark)
add_executable(benchhamming hamming_bench.cc)
target_link_libraries(benchhamming hamming benchmark::benchmark)

# Build tests
include(CTest)
if(BUILD_TESTING)
  add_subdirectory(Catch2)
  include(Catch2/contrib/Catch.cmake)
  add_executable(test_hamming hamming_t.cc)
  target_link_libraries(test_hamming hamming Catch2::Catch2)
  catch_discover_tests(test_hamming)
endif()

# The standalone executable
add_executable(distance distance.cc)
target_link_libraries(distance PUBLIC hamming)

# The python package
add_subdirectory(pybind11)
pybind11_add_module(hammingdist EXCLUDE_FROM_ALL pythonlib.cc)
target_link_libraries(hammingdist PUBLIC hamming)

include(FeatureSummary)
feature_summary(WHAT ALL)
