#include "hamming/hamming.hh"
#include "hamming_impl.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace hamming {

DataSet::DataSet(std::vector<std::string> &data, bool include_x,
                 bool clear_input_data, std::vector<std::size_t> &&indices)
    : nsamples(data.size()), sequence_indices(std::move(indices)) {
  validate_data(data);
  result = distances(data, include_x, clear_input_data);
}

DataSet::DataSet(const std::string &filename) {
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

static std::size_t uint_sqrt(std::size_t x) {
  return static_cast<std::size_t>(std::round(std::sqrt(x)));
}

DataSet::DataSet(std::vector<DistIntType> &&distances)
    : result{std::move(distances)} {
  // infer n from number of lower triangular matrix elements = n(n-1)/2
  nsamples = (uint_sqrt(8 * result.size() + 1) + 1) / 2;
}

void DataSet::dump(const std::string &filename) {
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

void DataSet::dump_lower_triangular(const std::string &filename) {
  std::ofstream stream(filename);
#ifdef HAMMING_WITH_OPENMP
  std::size_t block_size = 200;
#pragma omp parallel for ordered schedule(static, 1)
  for (std::size_t i_start = 1; i_start < nsamples; i_start += block_size) {
    std::stringstream line;
    std::size_t i_end = i_start + block_size;
    if (i_end > nsamples) {
      i_end = nsamples;
    }
    for (std::size_t i = i_start; i < i_end; ++i) {
      std::size_t offset{i * (i - 1) / 2};
      for (std::size_t j = 0; j + 1 < i; ++j) {
        line << static_cast<int>(result[offset + j]) << ",";
      }
      line << static_cast<int>(result[offset + i - 1]) << "\n";
    }
#pragma omp ordered
    stream << line.str();
  }
#else
  std::size_t k = 0;
  for (std::size_t i = 1; i < nsamples; ++i) {
    for (std::size_t j = 0; j + 1 < i; ++j) {
      stream << static_cast<int>(result[k++]) << ",";
    }
    stream << static_cast<int>(result[k++]) << "\n";
  }
#endif
}

void DataSet::dump_sequence_indices(const std::string &filename) {
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

int DataSet::operator[](const std::array<std::size_t, 2> &index) const {
  auto i = index[0];
  auto j = index[1];
  if (i < j)
    return result[j * (j - 1) / 2 + i];
  if (i > j)
    return result[i * (i - 1) / 2 + j];
  // This is a diagonal entry
  return 0;
}

DataSet from_stringlist(std::vector<std::string> &data) {
  return DataSet(data);
}

DataSet from_csv(const std::string &filename) { return DataSet(filename); }

DataSet from_lower_triangular(const std::string &filename) {
  std::vector<DistIntType> distances;
  std::ifstream stream(filename);
  std::string line;
  while (std::getline(stream, line)) {
    std::istringstream s(line);
    std::string d;
    while (s.good()) {
      std::getline(s, d, ',');
      distances.push_back(safe_int_cast(std::stoi(d)));
    }
  }
  return DataSet(std::move(distances));
}

DataSet from_fasta(const std::string &filename, bool include_x,
                   bool remove_duplicates, std::size_t n) {
  std::vector<std::string> data;
  data.reserve(n);
  if (n == 0) {
    n = std::numeric_limits<std::size_t>::max();
    data.reserve(65536);
  }
  std::unordered_map<std::string, std::size_t> map_seq_to_index;
  std::vector<std::size_t> sequence_indices{};
  // Initializing the stream
  std::ifstream stream(filename);
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
    if (remove_duplicates) {
      auto result = map_seq_to_index.emplace(std::move(seq), count_unique);
      if (result.second) {
        ++count_unique;
      }
      sequence_indices.push_back(result.first->second);
    } else {
      data.push_back(std::move(seq));
    }
    ++count;
  }
  if (remove_duplicates) {
    // copy each unique sequence to the vector of strings
    data.resize(count_unique);
    for (auto &key_value_pair : map_seq_to_index) {
      data[key_value_pair.second] = key_value_pair.first;
    }
  }
  return DataSet(data, include_x, true, std::move(sequence_indices));
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

} // namespace hamming
