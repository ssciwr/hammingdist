cmake_minimum_required(VERSION 3.11..3.27)

project(
  hammingdist
  VERSION 0.21.0
  LANGUAGES CXX)

include(CTest)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 17)

# compilation options
set(HAMMING_WITH_OPENMP
    no
    CACHE BOOL "Build with OpenMP support (multithreading)")

set(HAMMING_BUILD_BENCHMARKS
    yes
    CACHE BOOL "Build hamming benchmarks")

set(HAMMING_BUILD_PYTHON
    yes
    CACHE BOOL "Build hammingdist python interface")

set(HAMMING_WITH_SSE2
    yes
    CACHE BOOL "Enable SSE2 optimized code on x86_64 CPUs")

set(HAMMING_WITH_AVX2
    yes
    CACHE BOOL "Enable AVX2 optimized code on x86_64 CPUs")

set(HAMMING_WITH_AVX512
    yes
    CACHE BOOL "Enable AVX512 optimized code on x86_64 CPUs")

set(HAMMING_WITH_NEON
    no
    CACHE BOOL "Enable NEON optimized code on Arm64 CPUs")

# Add git submodules
add_subdirectory(ext)

# Build the hamming library
add_subdirectory(src)

# The standalone executable
add_executable(distance distance.cc)
target_link_libraries(distance PUBLIC hamming)

# The python package
if(HAMMING_BUILD_PYTHON)
  add_subdirectory(python)
endif()

include(FeatureSummary)
feature_summary(WHAT ALL)
