#include "hamming/distance_neon.hh"
#include <arm_neon.h>

namespace hamming {

int distance_neon(const std::vector<GeneBlock> &a,
                  const std::vector<GeneBlock> &b, int max_dist) {
  // distance implementation using NEON simd intrinsics
  // a 128-bit register holds 16 GeneBlocks, i.e. 32 genes
  constexpr std::size_t n_geneblocks{16};
  int r{0};
  // mask to select LSB of each gene
  const uint8x16_t lsb = vdupq_n_u8(1);
  // mask to select lower gene from each GeneBlock
  const uint8x16_t mask0 = vdupq_n_u8(mask_gene0);
  // mask to select upper gene from each GeneBlock
  const uint8x16_t mask1 = vdupq_n_u8(mask_gene1);
  // vector of partial distance counts
  uint8x16_t r_s;
  // work registers
  uint8x16_t r_a;
  uint8x16_t r_b;
  // each iteration processes 16 GeneBlocks
  std::size_t n_iter{a.size() / n_geneblocks};
  // each partial distance count is stored in a uint8, so max value = 255,
  // and the value can be increased by at most 2 with each iteration,
  // so up to 127 inner iterations for a max value of 254 avoid overflow.
  // if max_dist is large then we maximise the number of inner iterations, but
  // if it is small then we do fewer inner iterations to allow more
  // opportunities to return early if we have reached max_dist
  std::size_t n_inner = max_dist >= 255 ? 127 : 16;
  std::size_t n_outer{1 + n_iter / n_inner};
  for (std::size_t j = 0; j < n_outer; ++j) {
    std::size_t n{std::min((j + 1) * n_inner, n_iter)};
    r_s = vdupq_n_u8(0);
    for (std::size_t i = j * n_inner; i < n; ++i) {
      // load a[i], b[i] into registers
      r_a = vld1q_u8(a.data() + n_geneblocks * i);
      r_b = vld1q_u8(b.data() + n_geneblocks * i);
      // a[i] & b[i]
      r_a = vandq_u8(r_a, r_b);
      // mask lower genes
      r_b = vandq_u8(r_a, mask0);
      // mask upper genes
      r_a = vandq_u8(r_a, mask1);
      // compare genes with zero to get either 00000000 or 11111111
      r_a = vceqzq_u8(r_a);
      r_b = vceqzq_u8(r_b);
      // only keep LSB for each uint8 to get either 0 or 1
      r_a = vandq_u8(r_a, lsb);
      r_b = vandq_u8(r_b, lsb);
      // add these values to distance counts
      r_s = vaddq_u8(r_s, r_a);
      r_s = vaddq_u8(r_s, r_b);
    }
    // sum the 16 distances in r_s & add to r
    r += vaddlvq_u8(r_s);
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
