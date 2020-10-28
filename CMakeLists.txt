project(covid-tda LANGUAGES CXX)

cmake_minimum_required(VERSION 3.11)

# Build library
add_library(hamming STATIC hamming.cc)
target_compile_options(hamming PUBLIC -fPIC)

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
