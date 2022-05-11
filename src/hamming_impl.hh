#ifndef _HAMMING_IMPL_HH
#define _HAMMING_IMPL_HH

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "hamming/hamming_types.hh"
#include "hamming_impl_types.hh"

namespace hamming {

std::array<Gene, 256> lookupTable(bool include_x = false);

DistIntType safe_int_cast(int x);

std::vector<DistIntType> distances(std::vector<std::string> &data,
                                   bool include_x, bool clear_input_data);

int distance_sparse(const SparseData &a, const SparseData &b);

int distance_cpp(const std::vector<Gene> &a, const std::vector<Gene> &b);

void validate_data(const std::vector<std::string> &data);

std::vector<SparseData> to_sparse_data(const std::vector<std::string> &data,
                                       bool include_x);

std::vector<std::vector<Gene>>
to_dense_data(const std::vector<std::string> &data);

std::vector<Gene> from_string(const std::string &str);

} // namespace hamming

#endif
