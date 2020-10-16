project(covid-tda LANGUAGES CXX)

cmake_minimum_required(VERSION 3.11)

find_package(OpenMP)

add_executable(distance distance.cc)

add_subdirectory(alternatives)

include(FeatureSummary)
feature_summary(WHAT ALL)
