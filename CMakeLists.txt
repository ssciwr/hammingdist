cmake_minimum_required(VERSION 3.23..3.30)

project(
  hammingdist
  VERSION 1.3.0
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

set(HAMMING_WITH_CUDA
    no
    CACHE BOOL "Build with CUDA support (nvidia gpu)")

if(HAMMING_WITH_CUDA)
  enable_language(CUDA)
endif()

# Add git submodules
add_subdirectory(ext)

# Build the hamming library
add_subdirectory(src)

# The python package
if(HAMMING_BUILD_PYTHON)
  add_subdirectory(python)
endif()

include(FeatureSummary)
feature_summary(WHAT ALL)
