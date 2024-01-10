add_subdirectory(fmt)

if(HAMMING_BUILD_PYTHON)
  add_subdirectory(pybind11)
endif()

if(BUILD_TESTING)
  add_subdirectory(Catch2)
endif()

if(HAMMING_BUILD_BENCHMARKS)
  set(BENCHMARK_ENABLE_GTEST_TESTS OFF)
  set(BENCHMARK_ENABLE_TESTING OFF)
  add_subdirectory(benchmark)
endif()

# temporarily set BUILD_TESTING=OFF for cpu_features
set(OUR_BUILD_TESTING_VALUE "${BUILD_TESTING}")
set(BUILD_TESTING OFF)
set(BUILD_PIC ON)
add_subdirectory(cpu_features)
# restore our original value of BUILD_TESTING
set(BUILD_TESTING "${OUR_BUILD_TESTING_VALUE}")
unset(OUR_BUILD_TESTING_VALUE)
