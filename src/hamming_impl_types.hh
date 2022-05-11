#ifndef _HAMMING_IMPL_TYPES_HH
#define _HAMMING_IMPL_TYPES_HH

#include <array>
#include <cstdint>

namespace hamming {

// 8-bit representation of gene:
using Gene = std::uint8_t;
using SparseData = std::vector<std::size_t>;
constexpr std::size_t n_bits_per_gene{8};

} // namespace hamming

#endif
