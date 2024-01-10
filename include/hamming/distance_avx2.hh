#pragma once

#include <cstdint>
#include <vector>

#include "hamming/hamming_impl_types.hh"

namespace hamming {

int distance_avx2(const std::vector<GeneBlock> &a,
                  const std::vector<GeneBlock> &b);

}
