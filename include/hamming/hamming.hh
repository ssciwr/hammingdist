#pragma once

#include "hamming/hamming_impl.hh"
#include "hamming/hamming_types.hh"
#include "hamming/hamming_utils.hh"
#ifdef HAMMING_WITH_OPENMP
#include <omp.h>
#endif
#include <cmath>
#include <fmt/core.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace hamming {

inline std::size_t uint_sqrt(std::size_t x) {
  return static_cast<std::size_t>(std::round(std::sqrt(x)));
}

template <typename DistIntType> struct DataSet {
  explicit DataSet(std::vector<std::string> &data, bool include_x = false,
                   bool clear_input_data = false,
                   std::vector<std::size_t> &&indices = {},
                   bool use_gpu = false,
                   int max_distance = std::numeric_limits<int>::max())
      : nsamples(data.size()), sequence_indices(std::move(indices)) {
    validate_data(data);
    result = distances<DistIntType>(data, include_x, clear_input_data, use_gpu,
                                    max_distance);
  }

  explicit DataSet(const std::string &filename) {
    // Determine correct dataset size
    std::ifstream stream(filename);
    std::string line;
    nsamples = std::count(std::istreambuf_iterator<char>(stream),
                          std::istreambuf_iterator<char>(), '\n');
    result.resize(nsamples * (nsamples + 1) / 2);

    // Read the data
    stream = std::ifstream(filename);
    std::size_t i{0};
    std::size_t current{0};
    while (std::getline(stream, line)) {
      std::istringstream s(line);
      std::string d;
      for (std::size_t j = 0; j < current; ++j) {
        std::getline(s, d, ',');
        result[i++] = std::stoi(d);
      }
      ++current;
    }
  }

  explicit DataSet(std::vector<DistIntType> &&distances)
      : result{std::move(distances)} {
    // infer n from number of lower triangular matrix elements = n(n-1)/2
    nsamples = (uint_sqrt(8 * result.size() + 1) + 1) / 2;
  }

  void dump(const std::string &filename) {
    std::ofstream stream(filename);
    for (std::size_t i = 0; i < nsamples; ++i) {
      for (std::size_t j = 0; j < nsamples; ++j) {
        stream << (*this)[{i, j}];
        if (j != nsamples - 1)
          stream << ", ";
      }
      stream << std::endl;
    }
  }

  void dump_lower_triangular(const std::string &filename) {
    write_lower_triangular(filename, result);
  }

  void dump_sparse(const std::string &filename, int threshold) {
    std::ofstream stream(filename);
#ifdef HAMMING_WITH_OPENMP
    constexpr std::size_t samples_per_thread{200};
#pragma omp parallel for schedule(static, 1) default(none)                     \
    shared(result, stream, nsamples, samples_per_thread, threshold)
    for (std::size_t i_start = 1; i_start < nsamples;
         i_start += samples_per_thread) {
      std::size_t i_end{std::min(i_start + samples_per_thread, nsamples)};
      std::size_t offset{i_start * (i_start - 1) / 2};
      auto *d = result.data() + offset;
      std::string thread_local_lines;
      for (std::size_t i = i_start; i < i_end; ++i) {
        for (std::size_t j = 0; j < i; ++j) {
          if (*d <= threshold) {
            thread_local_lines.append(
                fmt::format("{} {} {}\n", i, j, static_cast<int>(*d)));
          }
          ++d;
        }
      }
#pragma omp critical
      stream << thread_local_lines;
    }
#else
    auto *d = result.data();
    for (std::size_t i = 1; i < nsamples; ++i) {
      for (std::size_t j = 0; j < i; ++j) {
        if (*d <= threshold) {
          stream << fmt::format("{} {} {}\n", i, j, static_cast<int>(*d));
        }
        ++d;
      }
    }
#endif
  }

  void dump_sequence_indices(const std::string &filename) {
    std::ofstream stream(filename);
    if (sequence_indices.empty()) {
      for (std::size_t i = 0; i < nsamples; ++i) {
        stream << i << "\n";
      }
      return;
    }
    for (auto sequence_index : sequence_indices) {
      stream << sequence_index << "\n";
    }
  }

  int operator[](const std::array<std::size_t, 2> &index) const {
    auto i = index[0];
    auto j = index[1];
    if (i < j)
      return result[j * (j - 1) / 2 + i];
    if (i > j)
      return result[i * (i - 1) / 2 + j];
    // This is a diagonal entry
    return 0;
  }

  std::size_t nsamples;
  std::vector<DistIntType> result;
  std::vector<std::size_t> sequence_indices{};
};

DataSet<DefaultDistIntType>
from_stringlist(std::vector<std::string> &data, bool include_x = false,
                bool use_gpu = false,
                int max_distance = std::numeric_limits<int>::max());

DataSet<DefaultDistIntType> from_csv(const std::string &filename);

template <typename DistIntType>
DataSet<DistIntType> from_lower_triangular(const std::string &filename) {
  std::vector<DistIntType> distances;
  std::ifstream stream(filename);
  std::string line;
  while (std::getline(stream, line)) {
    std::istringstream s(line);
    std::string d;
    while (s.good()) {
      std::getline(s, d, ',');
      distances.push_back(safe_int_cast<DistIntType>(std::stoi(d)));
    }
  }
  return DataSet(std::move(distances));
}

template <typename DistIntType>
DataSet<DistIntType>
from_fasta(const std::string &filename, bool include_x = false,
           bool remove_duplicates = false, std::size_t n = 0,
           bool use_gpu = false,
           int max_distance = std::numeric_limits<int>::max()) {
  auto [data, sequence_indices] = read_fasta(filename, remove_duplicates, n);
  return DataSet<DistIntType>(data, include_x, true,
                              std::move(sequence_indices), use_gpu,
                              max_distance);
}

void from_fasta_to_lower_triangular(
    const std::string &input_filename, const std::string &output_filename,
    bool remove_duplicates = false, std::size_t n = 0, bool use_gpu = false,
    int max_distance = std::numeric_limits<int>::max());

ReferenceDistIntType distance(const std::string &seq0, const std::string &seq1,
                              bool include_x = false);

std::vector<ReferenceDistIntType>
fasta_reference_distances(const std::string &reference_sequence,
                          const std::string &fasta_file,
                          bool include_x = false);

std::vector<std::size_t> fasta_sequence_indices(const std::string &fasta_file,
                                                std::size_t n = 0);

} // namespace hamming
