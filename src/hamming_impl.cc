#include "hamming/hamming_impl.hh"
#include <algorithm>
#if !(defined(__aarch64__) || defined(_M_ARM64))
#include <cpuinfo_x86.h>
#endif
#include <stdexcept>
#include <unordered_map>
#ifdef HAMMING_WITH_SSE2
#include "hamming/distance_sse2.hh"
#endif
#ifdef HAMMING_WITH_AVX2
#include "hamming/distance_avx2.hh"
#endif
#ifdef HAMMING_WITH_AVX512
#include "hamming/distance_avx512.hh"
#endif
#ifdef HAMMING_WITH_NEON
#include "hamming/distance_neon.hh"
#endif

namespace hamming {

// bit meaning:
// 1111: '-'
// 0001: 'A'
// 0010: 'C'
// 0100: 'G'
// 1000: 'T'
// 0000: invalid
// 0011: 'X' (only if include_x = true)
std::array<GeneBlock, 256> lookupTable(bool include_x) {
  std::array<GeneBlock, 256> lookup{};
  lookup[std::size_t('-')] = 0xff;
  lookup[std::size_t('A')] = 1 | (1 << n_bits_per_gene);
  lookup[std::size_t('C')] = (1 << 1) | (1 << (n_bits_per_gene + 1));
  lookup[std::size_t('G')] = (1 << 2) | (1 << (n_bits_per_gene + 2));
  lookup[std::size_t('T')] = (1 << 3) | (1 << (n_bits_per_gene + 3));
  if (include_x) {
    lookup[std::size_t('X')] =
        1 | (1 << 1) | (1 << n_bits_per_gene) | (1 << (n_bits_per_gene + 1));
  }

  return lookup;
}

distance_func_ptr get_fastest_supported_distance_func() {
  std::string simd_str = "no";
  distance_func_ptr distance_func{distance_cpp};
#if defined(__aarch64__) || defined(_M_ARM64)
#ifdef HAMMING_WITH_NEON
  distance_func = distance_neon;
  simd_str = "NEON";
#endif
#else
  const auto features = cpu_features::GetX86Info().features;
#ifdef HAMMING_WITH_SSE2
  if (features.sse2) {
    distance_func = distance_sse2;
    simd_str = "SSE2";
  }
#endif
#ifdef HAMMING_WITH_AVX2
  if (features.avx2) {
    distance_func = distance_avx2;
    simd_str = "AVX2";
  }
#endif
#ifdef HAMMING_WITH_AVX512
  if (features.avx512bw) {
    distance_func = distance_avx512;
    simd_str = "AVX512";
  }
#endif
#endif
  std::cout << "# hammingdist :: Using CPU with " << simd_str
            << " SIMD extensions..." << std::endl;
  return distance_func;
}

void validate_data(const std::vector<std::string> &data) {
  if (data.empty() || data[0].empty()) {
    throw std::runtime_error("Error: Empty sequence");
  }
  auto length{data[0].size()};
  for (const auto &d : data) {
    if (d.size() != length) {
      throw std::runtime_error(
          "Error: Sequences do not all have the same length");
    }
  }
}

int distance_sparse(const SparseData &a, const SparseData &b, int max_dist) {
  int r{0};
  std::size_t ia{0};
  std::size_t ib{0};
  while (ia < a.size() && ib < b.size()) {
    if (a[ia] < b[ib]) {
      // add distance contribution from a
      r += static_cast<int>(a[ia + 1] != 0xff);
      ia += 2;
    } else if (a[ia] > b[ib]) {
      // add distance contribution from b
      r += static_cast<int>(b[ib + 1] != 0xff);
      ib += 2;
    } else {
      // a and b have the same index: i.e. overlap
      // so add distance contribution from combination
      bool nodash = (a[ia + 1] != 0xff) && (b[ib + 1] != 0xff);
      bool differ = (a[ia + 1] != b[ib + 1]);
      r += static_cast<int>(nodash && differ);
      // advance both
      ia += 2;
      ib += 2;
    }
    if (r >= max_dist) {
      return max_dist;
    }
  }
  while (ia < a.size()) {
    r += static_cast<int>(a[ia + 1] != 0xff);
    ia += 2;
  }
  while (ib < b.size()) {
    r += static_cast<int>(b[ib + 1] != 0xff);
    ib += 2;
  }
  return std::min(r, max_dist);
}

int distance_cpp(const std::vector<GeneBlock> &a,
                 const std::vector<GeneBlock> &b, int max_dist) {
  int r{0};
  for (std::size_t i = 0; i < a.size(); ++i) {
    auto c{static_cast<GeneBlock>(a[i] & b[i])};
    r += static_cast<int>((c & mask_gene0) == 0) +
         static_cast<int>((c & mask_gene1) == 0);
  }
  return std::min(r, max_dist);
}

static std::string
get_reference_expression(const std::vector<std::string> &data, bool include_x) {
  std::string g0;
  std::size_t length{data[0].size()};
  g0.reserve(length);
  std::array<std::size_t, 256> ctoi{0};
  ctoi[static_cast<std::size_t>('A')] = 1;
  ctoi[static_cast<std::size_t>('C')] = 2;
  ctoi[static_cast<std::size_t>('G')] = 3;
  ctoi[static_cast<std::size_t>('T')] = 4;
  if (include_x) {
    ctoi[static_cast<std::size_t>('X')] = 5;
  }
  std::vector<std::array<std::size_t, 6>> counts(length,
                                                 std::array<std::size_t, 6>{});
  for (const auto &g : data) {
    for (std::size_t i = 0; i < g.size(); ++i) {
      ++(counts[i][ctoi[g[i]]]);
    }
  }
  std::array<char, 6> itoc{'A', 'A', 'C', 'G', 'T', 'X'};
  for (const auto &count : counts) {
    g0.push_back(itoc[std::distance(
        count.cbegin(), std::max_element(count.cbegin() + 1, count.cend()))]);
  }
  return g0;
}

std::vector<SparseData> to_sparse_data(const std::vector<std::string> &data,
                                       bool include_x) {
  std::vector<SparseData> sparseData;
  sparseData.reserve(data.size());
  auto lookup = lookupTable(include_x);
  auto seq0 = get_reference_expression(data, include_x);
  std::size_t count{0};
  for (const auto &seq : data) {
    sparseData.emplace_back();
    auto &d = sparseData.back();
    for (std::size_t i = 0; i < seq.size(); ++i) {
      if (seq0[i] != seq[i]) {
        d.push_back(i);
        d.push_back(lookup[seq[i]]);
      }
    }
  }
  return sparseData;
}

std::vector<std::vector<GeneBlock>>
to_dense_data(const std::vector<std::string> &data) {
  std::vector<std::vector<GeneBlock>> dense;
  dense.reserve(data.size());
  for (const auto &str : data) {
    dense.push_back(from_string(str));
  }
  return dense;
}

std::pair<std::vector<std::string>, std::vector<std::size_t>>
read_fasta(const std::string &filename, bool remove_duplicates, std::size_t n) {
  std::pair<std::vector<std::string>, std::vector<std::size_t>>
      data_and_sequence_indices;
  auto &[data, sequence_indices] = data_and_sequence_indices;
  data.reserve(n);
  if (n == 0) {
    n = std::numeric_limits<std::size_t>::max();
    data.reserve(65536);
  }
  std::unordered_map<std::string, std::size_t> map_seq_to_index;
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
  return data_and_sequence_indices;
}

std::vector<GeneBlock> from_string(const std::string &str) {
  alignas(16) std::vector<GeneBlock> r;
  auto lookup = lookupTable();
  std::size_t n_full_blocks{str.size() / 2};
  r.reserve(1 + n_full_blocks + 3);
  auto iter_str = str.cbegin();
  for (std::size_t i_block = 0; i_block < n_full_blocks; ++i_block) {
    r.push_back(lookup[*iter_str] & mask_gene0);
    ++iter_str;
    r.back() |= (lookup[*iter_str] & mask_gene1);
    ++iter_str;
  }
  // pad last GeneBlock if odd number of chars
  if (iter_str != str.cend()) {
    r.push_back(lookup[*iter_str] & mask_gene0);
    r.back() |= (lookup['-'] & mask_gene1);
  }
  // pad to ensure 64-bit alignment
  while (8 * (r.size() / 8) != r.size()) {
    r.push_back(lookup['-']);
  }

  return r;
}

} // namespace hamming
