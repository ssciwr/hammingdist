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

std::vector<int> distances(const std::vector<std::vector<GeneBlock>>& data){
  std::vector<int> result((data.size() - 1) * data.size()/2, 0);
  // choose fastest supported distance function
  const auto features = cpu_features::GetX86Info().features;
  int (*distance_func)(const std::vector<GeneBlock>& a, const std::vector<GeneBlock>& b) = distance_cpp;
#ifdef HAMMING_WITH_AVX512
  if(features.avx512bw){
    distance_func = distance_avx512;
  }
#endif
#ifdef HAMMING_WITH_AVX2
  if(features.avx2){
    distance_func = distance_avx2;
  }
#endif
#ifdef HAMMING_WITH_SSE2
  if(features.sse2){
    distance_func = distance_sse2;
  }
#endif
  std::size_t nsamples{data.size()};
#ifdef HAMMING_WITH_OPENMP
  #pragma omp parallel for
#endif
  for(std::size_t i=0; i<nsamples; ++i){
    std::size_t offset{i * (i - 1) / 2};
    for(std::size_t j=0; j<i; ++j)
      result[offset + j] = distance_func(data[i], data[j]);
  }
  return result;
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
