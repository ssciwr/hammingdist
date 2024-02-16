#pragma once

#include <cstdint>
#include <vector>

#include "hamming/hamming_impl_types.hh"
#include <limits>

namespace hamming {

int distance_neon(const std::vector<GeneBlock> &a,
                  const std::vector<GeneBlock> &b,
                  int max_dist = std::numeric_limits<int>::max());

}
