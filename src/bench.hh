#pragma once

#include <benchmark/benchmark.h>
#include <random>
#include <string>
#include <vector>

namespace hamming {

std::string make_string(int64_t n, std::mt19937 &gen, bool include_dash = true);

void randomize_n(std::string &str, std::size_t n, std::mt19937 &gen);

std::vector<std::string> make_stringlist(int64_t n, int64_t n_samples,
                                         std::mt19937 &gen);

void write_fasta(const std::string &filename, const std::string &seq,
                 std::size_t n_seq, std::mt19937 &gen,
                 std::size_t randomise_every_n = 200);

template <typename DistIntType>
std::vector<DistIntType> make_distances(int64_t n, std::mt19937 &gen) {
  std::vector<DistIntType> v{};
  auto n_elements{n * (n - 1) / 2};
  v.reserve(n_elements);
  std::uniform_int_distribution<> distrib(
      0, std::numeric_limits<DistIntType>::max());
  for (int64_t i = 0; i < n_elements; ++i) {
    v.push_back(distrib(gen));
  }
  return v;
}

const char *const benchmark_tmp_output_file{"tmp.bench.output"};
const char *const benchmark_tmp_input_file{"tmp.bench.input"};

} // namespace hamming
