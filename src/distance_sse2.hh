#ifndef _HAMMING_DISTANCE_SSE2_HH
#define _HAMMING_DISTANCE_SSE2_HH

#include <cstdint>
#include <vector>

namespace hamming {

int distance_sse2(const std::vector<std::uint8_t>& a, const std::vector<std::uint8_t>& b);

}

#endif
