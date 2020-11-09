#include"hamming_impl.hh"

#include<limits>
#include<stdexcept>
#include<immintrin.h>
#ifdef HAMMING_WITH_OPENMP
#include<omp.h>
#endif

// bit meaning:
// 1111: '-'
// 0001: 'A'
// 0010: 'C'
// 0100: 'G'
// 1000: 'T'
// 0000: invalid
std::array<GeneBlock, 256> lookupTable()
{
  std::array<GeneBlock, 256> lookup{};
  lookup[std::size_t('-')] = 0xff;
  lookup[std::size_t('A')] = 1 | (1 << n_bits_per_gene);
  lookup[std::size_t('C')] = (1 << 1) | (1 << (n_bits_per_gene+1));
  lookup[std::size_t('G')] = (1 << 2) | (1 << (n_bits_per_gene+2));
  lookup[std::size_t('T')] = (1 << 3) | (1 << (n_bits_per_gene+3));

  return lookup;
}

int distance(const std::vector<GeneBlock>& a, const std::vector<GeneBlock>& b){
    int r{0};
    for (std::size_t i=0; i<a.size(); ++i) {
      auto c{static_cast<GeneBlock>(a[i] & b[i])};
      r += static_cast<int>((c & mask_gene0) == 0);
      r += static_cast<int>((c & mask_gene1) == 0);
    }
    return r;
}

int distance_sse2(const std::vector<GeneBlock>& a, const std::vector<GeneBlock>& b){
  // distance implementation using SSE2 simd intrinsics
  // a 128-bit register holds 16 GeneBlocks, i.e. 32 genes
  int r{0};
  // mask to select LSB of each gene
  const __m128i lsb = _mm_set1_epi8(1);
  // zero value for each gene to compare to
  const __m128i zeros = _mm_set1_epi8(0);
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
  std::size_t n_iter{a.size()/16};
  // each partial distance count is stored in a unit8, so max value = 255,
  // and the value can be increased by at most 2 with each iteration,
  // so we do 127 inner iterations for a max value of 254 to avoid overflow
  std::size_t n_inner{127};
  std::size_t n_outer{1 + n_iter/n_inner};
  for (std::size_t j=0; j<n_outer; ++j) {
    std::size_t n{std::min((j+1)*n_inner, n_iter)};
    r_s = _mm_set1_epi8(0);
    for (std::size_t i=j*n_inner; i<n; ++i) {
      // load a[i], b[i] into registers
      r_a = _mm_load_si128((__m128i*)(a.data()+16*i));
      r_b = _mm_load_si128((__m128i*)(b.data()+16*i));
      // a[i] & b[i]
      r_a = _mm_and_si128(r_a, r_b);
      // mask lower genes
      r_b = _mm_and_si128(r_a, mask0);
      // mask upper genes
      r_a = _mm_and_si128(r_a, mask1);
      // compare lower genes with zero
      r_b = _mm_cmpeq_epi8(r_b, zeros);
      // convert comparison to 1 if zero, 0 otherwise
      r_b = _mm_and_si128(r_b, lsb);
      // add this value to distance counts
      r_s = _mm_add_epi8(r_b, r_s);
      // repeat for upper genes
      r_a = _mm_cmpeq_epi8(r_a, zeros);
      r_a = _mm_and_si128(r_a, lsb);
      r_s = _mm_add_epi8(r_a, r_s);
    }
    // sum the 16 distances in r_s & add to r
    alignas(16) GeneBlock r_partial[16];
    _mm_store_si128((__m128i*)r_partial, r_s);
    for(std::size_t i=0; i<16; ++i){
        r += r_partial[i];
    }
  }
  // do last partial block without simd intrinsics
  for (std::size_t i=16*(a.size()/16); i<a.size(); ++i) {
    auto c{static_cast<GeneBlock>(a[i] & b[i])};
    r += static_cast<int>((c & mask_gene0) == 0);
    r += static_cast<int>((c & mask_gene1) == 0);
  }
  return r;
}

void validate_data(const std::vector<std::vector<GeneBlock>>& data){
    if(data.empty() || data[0].empty()){
        throw std::runtime_error("Error: Empty sequence");
    }
    auto length{data[0].size()};
    for(const auto& d : data){
        if(d.size() != length){
            throw std::runtime_error("Error: Sequences do not all have the same length");
        }
    }
}

std::vector<GeneBlock> from_string(const std::string& str){
  alignas(16) std::vector<GeneBlock> r;
  auto lookup = lookupTable();
  std::size_t n_full_blocks{str.size()/2};
  r.reserve(1+n_full_blocks);
  auto iter_str = str.cbegin();
  for (std::size_t i_block = 0; i_block < n_full_blocks; ++i_block) {
    r.push_back(lookup[*iter_str] & mask_gene0);
    ++iter_str;
    r.back() |= (lookup[*iter_str] & mask_gene1);
    ++iter_str;
  }
  // pad last GeneBlock if odd number of chars
  if(iter_str != str.cend()){
    r.push_back(lookup[*iter_str] & mask_gene0);
    r.back() |= (lookup['-'] & mask_gene1);
  }
  return r;
}
