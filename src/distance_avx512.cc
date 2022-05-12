#include "distance_avx512.hh"
#include <immintrin.h>

namespace hamming {

int distance_avx512(const std::vector<Gene> &a, const std::vector<Gene> &b) {
  // distance implementation using AVX512 simd intrinsics
  // a 512-bit register holds 64 Genes
  constexpr std::size_t n_genes{64};
  int r{0};
  // mask to select LSB of each gene
  const __m512i lsb = _mm512_set1_epi8(1);
  // vector of distance counts
  __m512i r_s;
  // work registers
  __m512i r_a;
  __m512i r_b;
  // mask register
  __mmask64 r_m;
  // each iteration processes 64 GeneBlocks
  std::size_t n_iter{a.size() / n_genes};
  // each partial distance count is stored in a uint8, so max value = 255,
  // and the value can be increased by at most 1 with each iteration,
  // so we do 255 inner iterations to avoid overflow
  std::size_t n_inner{255};
  std::size_t n_outer{1 + n_iter / n_inner};
  for (std::size_t j = 0; j < n_outer; ++j) {
    std::size_t n{std::min((j + 1) * n_inner, n_iter)};
    r_s = _mm512_set1_epi8(0);
    for (std::size_t i = j * n_inner; i < n; ++i) {
      // load a[i], b[i] into registers
      r_a = _mm512_loadu_si512((__m512i *)(a.data() + n_genes * i));
      r_b = _mm512_loadu_si512((__m512i *)(b.data() + n_genes * i));
      // if (a[i] & b[i]) != 0 put 1 in mask, otherwise 0 in mask
      r_m = _mm512_testn_epi8_mask(r_a, r_b);
      // add 1 for each mask bit to distance counts
      r_s = _mm512_mask_add_epi8(r_s, r_m, lsb, r_s);
    }
    // sum the 64 distances in r_s & add to r
    // sum 8 blocks of 8 uint8 distances to 8 uint64
    constexpr std::size_t n_partialsums{n_genes / 8};
    r_s = _mm512_sad_epu8(r_s, _mm512_setzero_si512());
    r += _mm512_reduce_add_epi64(r_s);
    //    alignas(64) std::uint64_t r_partial[n_partialsums];
    //    _mm512_storeu_si512((__m512i *)r_partial, r_s);
    //    for (std::size_t i = 0; i < n_partialsums; ++i) {
    //      r += r_partial[i];
    //    }
  }
  // do last partial block without simd intrinsics
  for (std::size_t i = n_genes * n_iter; i < a.size(); ++i) {
    r += static_cast<int>((a[i] & b[i]) == 0);
  }
  return r;
}

} // namespace hamming
