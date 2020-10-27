project(covid-tda LANGUAGES CXX)

cmake_minimum_required(VERSION 3.11)

# Build library
add_library(libhamming SHARED hamming.cc)

# The standalone executable
add_executable(distance distance.cc)
target_link_libraries(distance PUBLIC libhamming)

# The python package
add_subdirectory(pybind11)
pybind11_add_module(hammingdist pythonlib.cc)
target_link_libraries(hammingdist PUBLIC libhamming)

include(FeatureSummary)
feature_summary(WHAT ALL)
