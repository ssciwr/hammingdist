#ifndef _HAMMING_IMPL_HH
#define _HAMMING_IMPL_HH

#include<array>
#include<cstdint>
#include<string>
#include<vector>

namespace hamming {

// 4-bit representation of gene:
using GeneBlock = std::uint_fast8_t;
constexpr std::size_t n_bits_per_gene{4};
constexpr GeneBlock mask_gene0{0x0f};
constexpr GeneBlock mask_gene1{0xf0};

std::array<GeneBlock, 256> lookupTable();

std::vector<int> distances(const std::vector<std::vector<GeneBlock>>& data);

int distance_cpp(const std::vector<GeneBlock>& a, const std::vector<GeneBlock>& b);

int distance(const std::vector<GeneBlock>& a, const std::vector<GeneBlock>& b);

void validate_data(const std::vector<std::vector<GeneBlock>>& data);

std::vector<GeneBlock> from_string(const std::string& str);

}

#endif
