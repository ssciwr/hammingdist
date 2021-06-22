#ifndef _HAMMING_DISTANCE_SSE2_HH
#define _HAMMING_DISTANCE_SSE2_HH

#include <cstdint>
#include <vector>

#include "hamming_impl_types.hh"

namespace hamming {

int distance_sse2(const std::vector<GeneBlock>& a, const std::vector<GeneBlock>& b);

}

#endif
