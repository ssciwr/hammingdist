#ifndef _HAMMING_DISTANCE_AVX512_HH
#define _HAMMING_DISTANCE_AVX512_HH

#include <cstdint>
#include <vector>

#include "hamming_impl_types.hh"

namespace hamming {

int distance_avx512(const std::vector<GeneBlock> &a,
                    const std::vector<GeneBlock> &b);

}

#endif
