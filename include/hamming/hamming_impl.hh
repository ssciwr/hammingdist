#ifndef _HAMMING_IMPL_HH
#define _HAMMING_IMPL_HH

#include <array>
#if !(defined(__aarch64__) || defined(_M_ARM64))
#include <cpuinfo_x86.h>
#endif
#include <cstdint>
#include <limits>
#include <string>
#include <vector>
#ifdef HAMMING_WITH_OPENMP
#include <omp.h>
#endif
#ifdef HAMMING_WITH_NEON
#include "hamming/distance_neon.hh"
#endif
#ifdef HAMMING_WITH_SSE2
#include "hamming/distance_sse2.hh"
#endif
#ifdef HAMMING_WITH_AVX2
#include "hamming/distance_avx2.hh"
#endif
#ifdef HAMMING_WITH_AVX512
#include "hamming/distance_avx512.hh"
#endif

#include "hamming/hamming_impl_types.hh"
#include "hamming/hamming_types.hh"

namespace hamming {

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

std::vector<GeneBlock> from_string(const std::string &str);

template <typename DistIntType>
std::vector<DistIntType> distances(std::vector<std::string> &data,
                                   bool include_x, bool clear_input_data) {
  std::vector<DistIntType> result((data.size() - 1) * data.size() / 2, 0);
  auto sparse = to_sparse_data(data, include_x);
  std::size_t nsamples{data.size()};
  std::size_t sample_length{data[0].size()};

  // if X is included, we have to use the sparse distance function
  bool use_sparse = include_x;
  // otherwise, use heuristic to choose distance function: if < 0.5% of values
  // differ from reference genome, use sparse distance function
  if (!include_x) {
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
    if (clear_input_data) {
      data.clear();
    }
#ifdef HAMMING_WITH_OPENMP
#pragma omp parallel for
#endif
    for (std::size_t i = 0; i < nsamples; ++i) {
      std::size_t offset{i * (i - 1) / 2};
      for (std::size_t j = 0; j < i; ++j) {
        result[offset + j] =
            safe_int_cast<DistIntType>(distance_sparse(sparse[i], sparse[j]));
      }
    }
    return result;
  }

  // otherwise use fastest supported dense distance function
  auto dense = to_dense_data(data);
  if (clear_input_data) {
    data.clear();
  }
  int (*distance_func)(const std::vector<GeneBlock> &a,
                       const std::vector<GeneBlock> &b) = distance_cpp;

#if defined(__aarch64__) || defined(_M_ARM64)
#ifdef HAMMING_WITH_NEON
  distance_func = distance_neon;
#endif
#else
  const auto features = cpu_features::GetX86Info().features;
#ifdef HAMMING_WITH_SSE2
  if (features.sse2) {
    distance_func = distance_sse2;
  }
#endif
#ifdef HAMMING_WITH_AVX2
  if (features.avx2) {
    distance_func = distance_avx2;
  }
#endif
#ifdef HAMMING_WITH_AVX512
  if (features.avx512bw) {
    distance_func = distance_avx512;
  }
#endif
#endif

#ifdef HAMMING_WITH_OPENMP
#pragma omp parallel for schedule(static, 1)
#endif
  for (std::size_t i = 0; i < nsamples; ++i) {
    std::size_t offset{i * (i - 1) / 2};
    for (std::size_t j = 0; j < i; ++j)
      result[offset + j] =
          safe_int_cast<DistIntType>(distance_func(dense[i], dense[j]));
  }
  return result;
}

} // namespace hamming

#endif
