#ifndef _HAMMING_IMPL_TYPES_HH
#define _HAMMING_IMPL_TYPES_HH

#include<array>
#include<cstdint>

namespace hamming {

// 4-bit representation of gene:
using GeneBlock = std::uint8_t;
using SparseData = std::vector<std::size_t>;
constexpr std::size_t n_bits_per_gene{4};
constexpr GeneBlock mask_gene0{0x0f};
constexpr GeneBlock mask_gene1{0xf0};

}

#endif
