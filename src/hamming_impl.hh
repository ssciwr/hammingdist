#ifndef _HAMMING_IMPL_HH
#define _HAMMING_IMPL_HH

#include<array>
#include<cstdint>
#include<string>
#include<vector>

#include"hamming/hamming_types.hh"
#include"hamming_impl_types.hh"

namespace hamming {

std::array<GeneBlock, 256> lookupTable();

std::vector<DistIntType> distances(std::vector<std::string>& data, bool clear_input_data);

int distance_sparse(const SparseData& a, const SparseData& b);

int distance_cpp(const std::vector<GeneBlock>& a, const std::vector<GeneBlock>& b);

void validate_data(const std::vector<std::string>& data);

std::vector<SparseData> to_sparse_data(const std::vector<std::string>& data);

std::vector<std::vector<GeneBlock>> to_dense_data(const std::vector<std::string>& data);

std::vector<GeneBlock> from_string(const std::string& str);

}

#endif
