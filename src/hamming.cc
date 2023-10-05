#include "hamming/hamming.hh"
#include "hamming/hamming_impl.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

namespace hamming {

DataSet<DefaultDistIntType> from_stringlist(std::vector<std::string> &data,
                                            bool include_x, bool use_gpu) {
  return DataSet<DefaultDistIntType>(data, include_x, false, {}, use_gpu);
}

DataSet<DefaultDistIntType> from_csv(const std::string &filename) {
  return DataSet<DefaultDistIntType>(filename);
}

void from_fasta_to_lower_triangular(const std::string &input_filename,
                                    const std::string &output_filename,
                                    bool remove_duplicates, std::size_t n,
                                    bool use_gpu) {
  if (use_gpu) {
    std::cout << "# hammingdist :: Using GPU..." << std::endl;
  }
  auto start_time = std::chrono::high_resolution_clock::now();
  auto [data, sequence_indices] =
      read_fasta(input_filename, remove_duplicates, n);
  auto dense_data = to_dense_data(data);
  std::cout << "# hammingdist :: ...pre-processing completed in "
            << std::chrono::duration_cast<std::chrono::milliseconds>(
                   std::chrono::high_resolution_clock::now() - start_time)
                   .count()
            << " ms..." << std::endl;
#ifdef HAMMING_WITH_CUDA
  if (use_gpu) {
    distances_cuda_to_lower_triangular(dense_data, output_filename);
    return;
  }
#endif
  throw std::runtime_error(
      "from_fasta_to_lower_triangular is currently only available on GPU");
}

ReferenceDistIntType distance(const std::string &seq0, const std::string &seq1,
                              bool include_x) {
  auto lookup{lookupTable(include_x)};
  ReferenceDistIntType distance{0};
  if (seq0.size() != seq1.size()) {
    throw std::runtime_error(
        "Error: Sequences do not all have the same length");
  }
  for (std::size_t i = 0; i < seq0.size(); ++i) {
    auto a{lookup[seq0[i]]};
    auto b{lookup[seq1[i]]};
    bool invalid{(a & b) == 0x00};
    bool differ{(seq0[i] != seq1[i])};
    if (invalid || differ) {
      bool nodash{(a != 0xff) && (b != 0xff)};
      distance +=
          static_cast<ReferenceDistIntType>(invalid || (differ && nodash));
    }
  }
  return distance;
}

std::vector<std::size_t> fasta_sequence_indices(const std::string &fasta_file,
                                                std::size_t n) {
  std::vector<std::size_t> sequence_indices{};
  std::unordered_map<std::string, std::size_t> map_seq_to_index;
  std::ifstream stream(fasta_file);
  if (n == 0) {
    n = std::numeric_limits<std::size_t>::max();
  }
  std::size_t count = 0;
  std::size_t count_unique = 0;
  std::string line;
  // skip first header
  std::getline(stream, line);
  while (count < n && !stream.eof()) {
    std::string seq{};
    while (std::getline(stream, line) && line[0] != '>') {
      seq.append(line);
    }
    auto result = map_seq_to_index.emplace(std::move(seq), count_unique);
    if (result.second) {
      ++count_unique;
    }
    sequence_indices.push_back(result.first->second);
    ++count;
  }
  return sequence_indices;
}

std::vector<ReferenceDistIntType>
fasta_reference_distances(const std::string &reference_sequence,
                          const std::string &fasta_file, bool include_x) {
  std::vector<ReferenceDistIntType> distances;
  distances.reserve(65536);
  auto lookup{lookupTable(include_x)};
  std::vector<GeneBlock> ref;
  ref.reserve(reference_sequence.size());
  for (char c : reference_sequence) {
    ref.push_back(lookup[c]);
  }
  std::ifstream stream(fasta_file);
  if (!stream) {
    throw std::runtime_error("Error: Failed to open file '" + fasta_file + "'");
  }
  std::string line;
  // skip first header
  std::getline(stream, line);
  while (!stream.eof()) {
    std::string seq{};
    while (std::getline(stream, line) && line[0] != '>') {
      seq.append(line);
    }
    ReferenceDistIntType distance{0};
    for (std::size_t i = 0; i < seq.size(); ++i) {
      auto a{lookup[seq[i]]};
      bool invalid{(a & ref[i]) == 0x00};
      bool differ{(seq[i] != reference_sequence[i])};
      if (invalid || differ) {
        bool nodash{(a != 0xff) && (ref[i] != 0xff)};
        distance +=
            static_cast<ReferenceDistIntType>(invalid || (differ && nodash));
      }
    }
    distances.push_back(distance);
  }
  return distances;
}

} // namespace hamming
