#ifndef _HAMMING_DISTANCE_AVX512_HH
#define _HAMMING_DISTANCE_AVX512_HH

#include <cstdint>
#include <vector>

namespace hamming {

int distance_avx512(const std::vector<std::uint8_t>& a, const std::vector<std::uint8_t>& b);

}

#endif
