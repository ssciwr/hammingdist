#include "hamming/distance_sse2.hh"
#include <immintrin.h>

namespace hamming {

int distance_sse2(const std::vector<GeneBlock> &a,
                  const std::vector<GeneBlock> &b) {
  // distance implementation using SSE2 simd intrinsics
  // a 128-bit register holds 16 GeneBlocks, i.e. 32 genes
  constexpr std::size_t n_geneblocks{16};
  int r{0};
  // mask to select LSB of each gene
  const __m128i lsb = _mm_set1_epi8(1);
  // mask to select lower gene from each GeneBlock
  const __m128i mask0 = _mm_set1_epi8(mask_gene0);
  // mask to select upper gene from each GeneBlock
  const __m128i mask1 = _mm_set1_epi8(mask_gene1);
  // vector of distance counts
  __m128i r_s;
  // work registers
  __m128i r_a;
  __m128i r_b;
  // each iteration processes 16 GeneBlocks
  std::size_t n_iter{a.size() / n_geneblocks};
  // each partial distance count is stored in a unit8, so max value = 255,
  // and the value can be increased by at most 2 with each iteration,
  // so we do 127 inner iterations for a max value of 254 to avoid overflow
  std::size_t n_inner{127};
  std::size_t n_outer{1 + n_iter / n_inner};
  for (std::size_t j = 0; j < n_outer; ++j) {
    std::size_t n{std::min((j + 1) * n_inner, n_iter)};
    r_s = _mm_set1_epi8(0);
    for (std::size_t i = j * n_inner; i < n; ++i) {
      // load a[i], b[i] into registers
      r_a = _mm_load_si128((__m128i *)(a.data() + n_geneblocks * i));
      r_b = _mm_load_si128((__m128i *)(b.data() + n_geneblocks * i));
      // a[i] & b[i]
      r_a = _mm_and_si128(r_a, r_b);
      // mask lower genes
      r_b = _mm_and_si128(r_a, mask0);
      // mask upper genes
      r_a = _mm_and_si128(r_a, mask1);
      // compare lower genes with zero
      r_b = _mm_cmpeq_epi8(r_b, _mm_setzero_si128());
      // convert comparison to 1 if zero, 0 otherwise
      r_b = _mm_and_si128(r_b, lsb);
      // add this value to distance counts
      r_s = _mm_add_epi8(r_b, r_s);
      // repeat for upper genes
      r_a = _mm_cmpeq_epi8(r_a, _mm_setzero_si128());
      r_a = _mm_and_si128(r_a, lsb);
      r_s = _mm_add_epi8(r_a, r_s);
    }
    // sum the 16 distances in r_s & add to r
    // sum 2 blocks of 8 uint8 distances to 2 uint64
    constexpr std::size_t n_partialsums{n_geneblocks / 8};
    r_s = _mm_sad_epu8(r_s, _mm_setzero_si128());
    alignas(16) std::uint64_t r_partial[n_partialsums];
    _mm_store_si128((__m128i *)r_partial, r_s);
    for (std::size_t i = 0; i < n_partialsums; ++i) {
      r += r_partial[i];
    }
  }
  // do last partial block without simd intrinsics
  for (std::size_t i = n_geneblocks * n_iter; i < a.size(); ++i) {
    auto c{static_cast<GeneBlock>(a[i] & b[i])};
    r += static_cast<int>((c & mask_gene0) == 0);
    r += static_cast<int>((c & mask_gene1) == 0);
  }
  return r;
}

} // namespace hamming
