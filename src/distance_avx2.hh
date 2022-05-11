#ifndef _HAMMING_DISTANCE_AVX2_HH
#define _HAMMING_DISTANCE_AVX2_HH

#include <cstdint>
#include <vector>

#include "hamming_impl_types.hh"

namespace hamming {

int distance_avx2(const std::vector<Gene> &a, const std::vector<Gene> &b);

}

#endif
