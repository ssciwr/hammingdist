project(covid-tda LANGUAGES CXX)

cmake_minimum_required(VERSION 3.11)

# The standalone executable
add_executable(distance distance.cc)

# The python package
add_subdirectory(pybind11)
pybind11_add_module(hammingdist pythonlib.cc)


add_subdirectory(alternatives)

include(FeatureSummary)
feature_summary(WHAT ALL)
