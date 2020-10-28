project(covid-tda LANGUAGES CXX)

cmake_minimum_required(VERSION 3.11)

# Build library
add_library(hamming STATIC hamming.cc)
target_compile_options(hamming PUBLIC -fPIC)

# Build library benchmarks
set(BENCHMARK_ENABLE_GTEST_TESTS OFF)
add_subdirectory(benchmark)
add_executable(benchhamming hamming_bench.cc)
target_link_libraries(benchhamming libhamming benchmark::benchmark)

# The standalone executable
add_executable(distance distance.cc)
target_link_libraries(distance PUBLIC hamming)

# The python package
add_subdirectory(pybind11)
pybind11_add_module(hammingdist EXCLUDE_FROM_ALL pythonlib.cc)
target_link_libraries(hammingdist PUBLIC hamming)

include(FeatureSummary)
feature_summary(WHAT ALL)
