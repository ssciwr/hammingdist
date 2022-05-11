#include "distance_avx2.hh"
#include "hamming_impl_types.hh"
#include <immintrin.h>

namespace hamming {

int distance_avx2(const std::vector<Gene> &a, const std::vector<Gene> &b) {
  // distance implementation using AVX2 simd intrinsics
  // a 256-bit register holds 32 Genes
  constexpr std::size_t n_genes{32};
  int r{0};
  // mask to select LSB of each gene
  const __m256i lsb = _mm256_set1_epi8(1);
  // vector of distance counts
  __m256i r_s;
  // work registers
  __m256i r_a;
  __m256i r_b;
  // each iteration processes 32 Genes
  std::size_t n_iter{a.size() / n_genes};
  // each partial distance count is stored in a uint8, so max value = 255,
  // and the value can be increased by at most 1 with each iteration,
  // so we do 255 inner iterations to avoid overflow
  std::size_t n_inner{255};
  std::size_t n_outer{1 + n_iter / n_inner};
  for (std::size_t j = 0; j < n_outer; ++j) {
    std::size_t n{std::min((j + 1) * n_inner, n_iter)};
    r_s = _mm256_set1_epi8(0);
    for (std::size_t i = j * n_inner; i < n; ++i) {
      // load a[i], b[i] into registers
      r_a = _mm256_loadu_si256((__m256i *)(a.data() + n_genes * i));
      r_b = _mm256_loadu_si256((__m256i *)(b.data() + n_genes * i));
      // a[i] & b[i]
      r_a = _mm256_and_si256(r_a, r_b);
      // compare with zero: sets FF if equal, 00 otherwise
      r_b = _mm256_cmpeq_epi8(r_a, _mm256_setzero_si256());
      // convert FF to 1
      r_b = _mm256_and_si256(r_b, lsb);
      // add these 1's or 0's to the distance counts
      r_s = _mm256_add_epi8(r_b, r_s);
    }
    // sum the 32 distances in r_s & add to r
    // sum 4 blocks of 8 uint8 distances to 4 uint64
    constexpr std::size_t n_partialsums{n_genes / 8};
    r_s = _mm256_sad_epu8(r_s, _mm256_setzero_si256());
    alignas(32) std::uint64_t r_partial[n_partialsums];
    _mm256_storeu_si256((__m256i *)r_partial, r_s);
    for (std::size_t i = 0; i < n_partialsums; ++i) {
      r += r_partial[i];
    }
  }
  // do last partial block without simd intrinsics
  for (std::size_t i = n_genes * n_iter; i < a.size(); ++i) {
    r += static_cast<int>((a[i] & b[i]) == 0);
  }
  return r;
}

} // namespace hamming
