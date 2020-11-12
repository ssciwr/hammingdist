#include"hamming_impl.hh"
#ifdef HAMMING_WITH_SSE2
#include"distance_sse2.hh"
#endif
#ifdef HAMMING_WITH_AVX2
#include"distance_avx2.hh"
#endif
#ifdef HAMMING_WITH_AVX512
#include"distance_avx512.hh"
#endif
#include<limits>
#include<stdexcept>
#include<cpuinfo_x86.h>
#ifdef HAMMING_WITH_OPENMP
#include<omp.h>
#endif

namespace hamming {

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

int distance_cpp(const std::vector<GeneBlock>& a, const std::vector<GeneBlock>& b){
    int r{0};
    for (std::size_t i=0; i<a.size(); ++i) {
      auto c{static_cast<GeneBlock>(a[i] & b[i])};
      r += static_cast<int>((c & mask_gene0) == 0);
      r += static_cast<int>((c & mask_gene1) == 0);
    }
    return r;
}

int distance(const std::vector<GeneBlock>& a, const std::vector<GeneBlock>& b){
  const auto features = cpu_features::GetX86Info().features;
#ifdef HAMMING_WITH_AVX512
  if(features.avx512bw){
    return distance_avx512(a,b);
  }
#endif
#ifdef HAMMING_WITH_AVX2
  if(features.avx2){
    return distance_avx2(a,b);
  }
#endif
#ifdef HAMMING_WITH_SSE2
  if(features.sse2){
    return distance_sse2(a,b);
  }
#endif
  return distance_cpp(a,b);
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

}
