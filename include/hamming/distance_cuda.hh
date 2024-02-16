#pragma once

#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

#include "hamming/hamming_impl_types.hh"

namespace hamming {

bool distance_cuda_have_device();

int distance_cuda(const std::vector<GeneBlock> &a,
                  const std::vector<GeneBlock> &b,
                  int max_dist = std::numeric_limits<int>::max());

// for now explicit function def for each choice of integer type
std::vector<uint8_t>
distances_cuda_8bit(const std::vector<std::vector<GeneBlock>> &data,
                    uint8_t max_dist = std::numeric_limits<uint8_t>::max());

std::vector<uint16_t>
distances_cuda_16bit(const std::vector<std::vector<GeneBlock>> &data,
                     uint16_t max_dist = std::numeric_limits<uint16_t>::max());

void distances_cuda_to_lower_triangular(
    const std::vector<std::vector<GeneBlock>> &data,
    const std::string &filename, int max_distance);

} // namespace hamming
