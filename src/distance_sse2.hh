#ifndef _HAMMING_DISTANCE_SSE2_HH
#define _HAMMING_DISTANCE_SSE2_HH

#include <cstdint>
#include <vector>

#include "hamming_impl_types.hh"

namespace hamming {

int distance_sse2(const std::vector<Gene> &a, const std::vector<Gene> &b);

}

#endif
