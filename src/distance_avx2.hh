#ifndef _HAMMING_DISTANCE_AVX2_HH
#define _HAMMING_DISTANCE_AVX2_HH

#include <cstdint>
#include <vector>

#include "hamming_impl_types.hh"

namespace hamming {

int distance_avx2(const std::vector<GeneBlock>& a, const std::vector<GeneBlock>& b);

}

#endif
