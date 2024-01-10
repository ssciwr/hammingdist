#pragma once

#include <array>
#include <charconv>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#ifdef HAMMING_WITH_OPENMP
#include <omp.h>
#endif
#ifdef HAMMING_WITH_CUDA
#include "hamming/distance_cuda.hh"
#endif
#include "hamming/hamming_impl_types.hh"
#include "hamming/hamming_types.hh"

namespace hamming {

inline bool cuda_gpu_available() {
#ifdef HAMMING_WITH_CUDA
  return distance_cuda_have_device();
#else
  return false;
#endif
}

typedef int (*distance_func_ptr)(const std::vector<GeneBlock> &,
                                 const std::vector<GeneBlock> &);

distance_func_ptr get_fastest_supported_distance_func();

std::array<GeneBlock, 256> lookupTable(bool include_x = false);

template <typename DistIntType> DistIntType safe_int_cast(int x) {
  if (x > std::numeric_limits<DistIntType>::max()) {
    return std::numeric_limits<DistIntType>::max();
  }
  return static_cast<DistIntType>(x);
}

void validate_data(const std::vector<std::string> &data);

int distance_sparse(const SparseData &a, const SparseData &b);

int distance_cpp(const std::vector<GeneBlock> &a,
                 const std::vector<GeneBlock> &b);

std::vector<SparseData> to_sparse_data(const std::vector<std::string> &data,
                                       bool include_x);

std::vector<std::vector<GeneBlock>>
to_dense_data(const std::vector<std::string> &data);

std::pair<std::vector<std::string>, std::vector<std::size_t>>
read_fasta(const std::string &filename, bool remove_duplicates = false,
           std::size_t n = 0);

std::vector<GeneBlock> from_string(const std::string &str);

template <typename DistIntType>
std::vector<DistIntType> distances(std::vector<std::string> &data,
                                   bool include_x, bool clear_input_data,
                                   bool use_gpu) {
  std::vector<DistIntType> result((data.size() - 1) * data.size() / 2, 0);
  auto start_time = std::chrono::high_resolution_clock::now();
  auto print_timing = [&start_time](const std::string &event,
                                    bool final = false) {
    std::cout << "# hammingdist :: ..." << event << " completed in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     std::chrono::high_resolution_clock::now() - start_time)
                     .count()
              << " ms.";
    if (!final) {
      std::cout << "..";
    }
    std::cout << std::endl;
    start_time = std::chrono::high_resolution_clock::now();
  };

#ifndef HAMMING_WITH_CUDA
  if (use_gpu) {
    throw std::runtime_error("hammingdist was not compiled with GPU support, "
                             "please set use_gpu=False");
  }
#endif
  auto sparse = to_sparse_data(data, include_x);
  std::size_t nsamples{data.size()};
  std::size_t sample_length{data[0].size()};

  if (include_x && use_gpu) {
    throw std::runtime_error("use_gpu=True cannot be used if include_x=True, "
                             "please set use_gpu=False");
  }

  // if X is included, we have to use the sparse distance function
  bool use_sparse = include_x;
  // otherwise, use heuristic to choose distance function: if < 0.5% of values
  // differ from reference genome, and we're using the CPU, use sparse distance
  // function
  if (!include_x && !use_gpu) {
    constexpr double sparse_threshold{0.005};
    std::size_t n_diff{0};
    for (const auto &s : sparse) {
      n_diff += s.size() / 2;
    }
    double frac_diff{static_cast<double>(n_diff) /
                     static_cast<double>(nsamples * sample_length)};
    use_sparse = frac_diff < sparse_threshold;
  }
  if (use_sparse) {
    std::cout << "# hammingdist :: Using CPU with sparse distance function..."
              << std::endl;
    if (clear_input_data) {
      data.clear();
    }
    print_timing("pre-processing");
#ifdef HAMMING_WITH_OPENMP
#pragma omp parallel for default(none) shared(result, sparse, nsamples)
#endif
    for (std::size_t i = 0; i < nsamples; ++i) {
      std::size_t offset{i * (i - 1) / 2};
      for (std::size_t j = 0; j < i; ++j) {
        result[offset + j] =
            safe_int_cast<DistIntType>(distance_sparse(sparse[i], sparse[j]));
      }
    }
    print_timing("distance calculation", true);
    return result;
  }

  // otherwise use the fastest supported dense distance function
  auto dense = to_dense_data(data);
  if (clear_input_data) {
    data.clear();
  }

#ifdef HAMMING_WITH_CUDA
  if (use_gpu) {
    std::cout << "# hammingdist :: Using GPU..." << std::endl;
    print_timing("pre-processing");
    if constexpr (sizeof(DistIntType) == 1) {
      return distances_cuda_8bit(dense);
    } else if constexpr (sizeof(DistIntType) == 2) {
      return distances_cuda_16bit(dense);
    } else {
      throw std::runtime_error("No GPU implementation available");
    }
  }
#endif

  auto distance_func{get_fastest_supported_distance_func()};
  print_timing("pre-processing");
#ifdef HAMMING_WITH_OPENMP
#pragma omp parallel for schedule(static, 1) default(none)                     \
    shared(result, dense, nsamples, distance_func)
#endif
  for (std::size_t i = 0; i < nsamples; ++i) {
    std::size_t offset{i * (i - 1) / 2};
    for (std::size_t j = 0; j < i; ++j) {
      result[offset + j] =
          safe_int_cast<DistIntType>(distance_func(dense[i], dense[j]));
    }
  }
  print_timing("distance calculation", true);
  return result;
}

} // namespace hamming
