#ifndef _HAMMING_IMPL_HH
#define _HAMMING_IMPL_HH

#include<array>
#include<cstdint>
#include<string>
#include<vector>

namespace hamming {

// 4-bit representation of gene:
using GeneBlock = std::uint_fast8_t;
using SparseData = std::vector<std::size_t>;
constexpr std::size_t n_bits_per_gene{4};
constexpr GeneBlock mask_gene0{0x0f};
constexpr GeneBlock mask_gene1{0xf0};

std::array<GeneBlock, 256> lookupTable();

std::vector<uint16_t> distances(std::vector<std::string>& data, bool clear_input_data);

int distance_sparse(const SparseData& a, const SparseData& b);

int distance_cpp(const std::vector<GeneBlock>& a, const std::vector<GeneBlock>& b);

void validate_data(const std::vector<std::string>& data);

std::vector<SparseData> to_sparse_data(const std::vector<std::string>& data);

std::vector<std::vector<GeneBlock>> to_dense_data(const std::vector<std::string>& data);

std::vector<GeneBlock> from_string(const std::string& str);

}

#endif
