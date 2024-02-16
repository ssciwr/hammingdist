#include "hamming/distance_avx2.hh"
#include <immintrin.h>

namespace hamming {

int distance_avx2(const std::vector<GeneBlock> &a,
                  const std::vector<GeneBlock> &b, int max_dist) {
  // distance implementation using AVX2 simd intrinsics
  // a 256-bit register holds 32 GeneBlocks, i.e. 64 genes
  constexpr std::size_t n_geneblocks{32};
  int r{0};
  // mask to select LSB of each gene
  const __m256i lsb = _mm256_set1_epi8(1);
  // mask to select lower gene from each GeneBlock
  const __m256i mask0 = _mm256_set1_epi8(mask_gene0);
  // mask to select upper gene from each GeneBlock
  const __m256i mask1 = _mm256_set1_epi8(mask_gene1);
  // vector of distance counts
  __m256i r_s;
  // work registers
  __m256i r_a;
  __m256i r_b;
  // each iteration processes 32 GeneBlocks
  std::size_t n_iter{a.size() / n_geneblocks};
  // each partial distance count is stored in a unit8, so max value = 255,
  // and the value can be increased by at most 2 with each iteration,
  // so up to 127 inner iterations for a max value of 254 avoid overflow.
  // if max_dist is large then we maximise the number of inner iterations, but
  // if it is small then we do fewer inner iterations to allow more
  // opportunities to return early if we have reached max_dist
  std::size_t n_inner = max_dist >= 255 ? 127 : 16;
  std::size_t n_outer{1 + n_iter / n_inner};
  for (std::size_t j = 0; j < n_outer; ++j) {
    std::size_t n{std::min((j + 1) * n_inner, n_iter)};
    r_s = _mm256_set1_epi8(0);
    for (std::size_t i = j * n_inner; i < n; ++i) {
      // load a[i], b[i] into registers
      r_a = _mm256_loadu_si256((__m256i *)(a.data() + n_geneblocks * i));
      r_b = _mm256_loadu_si256((__m256i *)(b.data() + n_geneblocks * i));
      // a[i] & b[i]
      r_a = _mm256_and_si256(r_a, r_b);
      // mask lower genes
      r_b = _mm256_and_si256(r_a, mask0);
      // mask upper genes
      r_a = _mm256_and_si256(r_a, mask1);
      // compare lower genes with zero
      r_b = _mm256_cmpeq_epi8(r_b, _mm256_setzero_si256());
      // convert comparison to 1 if zero, 0 otherwise
      r_b = _mm256_and_si256(r_b, lsb);
      // add this value to distance counts
      r_s = _mm256_add_epi8(r_b, r_s);
      // repeat for upper genes
      r_a = _mm256_cmpeq_epi8(r_a, _mm256_setzero_si256());
      r_a = _mm256_and_si256(r_a, lsb);
      r_s = _mm256_add_epi8(r_a, r_s);
    }
    // sum the 32 distances in r_s & add to r
    // sum 4 blocks of 8 uint8 distances to 4 uint64
    constexpr std::size_t n_partialsums{n_geneblocks / 8};
    r_s = _mm256_sad_epu8(r_s, _mm256_setzero_si256());
    alignas(32) std::uint64_t r_partial[n_partialsums];
    _mm256_storeu_si256((__m256i *)r_partial, r_s);
    for (std::size_t i = 0; i < n_partialsums; ++i) {
      r += r_partial[i];
    }
    if (r >= max_dist) {
      return max_dist;
    }
  }
  // do last partial block without simd intrinsics
  for (std::size_t i = n_geneblocks * n_iter; i < a.size(); ++i) {
    auto c{static_cast<GeneBlock>(a[i] & b[i])};
    r += static_cast<int>((c & mask_gene0) == 0);
    r += static_cast<int>((c & mask_gene1) == 0);
  }
  return std::min(max_dist, r);
}

} // namespace hamming
