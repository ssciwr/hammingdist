#pragma once

#include <cstdint>
#include <iostream>
#include <vector>

#include "hamming/hamming_impl_types.hh"

namespace hamming {

bool distance_cuda_have_device();

int distance_cuda(const std::vector<GeneBlock> &a,
                  const std::vector<GeneBlock> &b);

// for now explicit function def for each choice of integer type
std::vector<uint8_t>
distances_cuda_8bit(const std::vector<std::vector<GeneBlock>> &data);

std::vector<uint16_t>
distances_cuda_16bit(const std::vector<std::vector<GeneBlock>> &data);

void distances_cuda_to_lower_triangular(
    const std::vector<std::vector<GeneBlock>> &data,
    const std::string &filename);

} // namespace hamming
