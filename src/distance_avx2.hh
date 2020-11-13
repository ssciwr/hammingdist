#ifndef _HAMMING_DISTANCE_AVX2_HH
#define _HAMMING_DISTANCE_AVX2_HH

#include <cstdint>
#include <vector>

namespace hamming {

int distance_avx2(const std::vector<std::uint8_t>& a, const std::vector<std::uint8_t>& b);

}

#endif
