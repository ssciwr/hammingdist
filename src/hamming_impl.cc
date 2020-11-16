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
#include<algorithm>
#include<cpuinfo_x86.h>
#ifdef HAMMING_WITH_OPENMP
#include<omp.h>
#endif
#include <iostream>
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

std::vector<int> distances(const std::vector<std::string>& data){
  std::vector<int> result((data.size() - 1) * data.size()/2, 0);
  auto sparse = to_sparse_data(data);
  std::size_t nsamples{data.size()};

  // if < 0.5% of values differ from reference genome, use sparse distance function
  constexpr double sparse_threshold{0.005};
  std::size_t n_diff{0};
  for(const auto& s : sparse){
      n_diff += s.size()/2;
  }
  double frac_diff{static_cast<double>(n_diff) / static_cast<double>(data.size()*data[0].size())};
  if(frac_diff < sparse_threshold){
      std::cout << "!! " << frac_diff << std::endl;
    #ifdef HAMMING_WITH_OPENMP
      #pragma omp parallel for
    #endif
      for(std::size_t i=0; i<nsamples; ++i){
        std::size_t offset{i * (i - 1) / 2};
        for(std::size_t j=0; j<i; ++j)
          result[offset + j] = distance_sparse(sparse[i], sparse[j]);
      }
      return result;
  }

  // otherwise use fastest supported dense distance function
  auto dense = to_dense_data(data);
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
#ifdef HAMMING_WITH_OPENMP
  #pragma omp parallel for
#endif
  for(std::size_t i=0; i<nsamples; ++i){
    std::size_t offset{i * (i - 1) / 2};
    for(std::size_t j=0; j<i; ++j)
      result[offset + j] = distance_func(dense[i], dense[j]);
  }
  return result;
}

int distance_sparse(const SparseData& a, const SparseData& b){
    int r{0};
    std::size_t ia{0};
    std::size_t ib{0};
    while(ia < a.size() && ib < b.size()){
        if(a[ia] < b[ib]){
            // add distance contribution from a
            r += static_cast<int>(a[ia+1] != 0xff);
            ia += 2;
        } else if(a[ia] > b[ib]){
            // add distance contribution from b
            r += static_cast<int>(b[ib+1] != 0xff);
            ib += 2;
        } else {
            // a and b have the same index: i.e. overlap
            // so add distance contribution from combination
            r += static_cast<int>((a[ia+1] & b[ib+1]) == 0);
            // advance both
            ia += 2;
            ib += 2;
        }
    }
    while(ia < a.size()){
        r += static_cast<int>(a[ia+1] != 0xff);
        ia += 2;
    }
    while(ib < b.size()){
        r += static_cast<int>(b[ib+1] != 0xff);
        ib += 2;
    }
    return r;
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

void validate_data(const std::vector<std::string>& data){
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

static std::string get_reference_expression(const std::vector<std::string>& data){
    std::string g0;
    std::size_t length{data[0].size()};
    g0.reserve(length);
    std::array<std::size_t, 256> ctoi{0};
    ctoi[static_cast<std::size_t>('A')] = 1;
    ctoi[static_cast<std::size_t>('C')] = 2;
    ctoi[static_cast<std::size_t>('G')] = 3;
    ctoi[static_cast<std::size_t>('T')] = 4;
    std::vector<std::array<std::size_t, 5>> counts(length, std::array<std::size_t, 5>{});
    for (const auto& g : data){
        for (std::size_t i=0; i<g.size(); ++i){
            ++(counts[i][ctoi[g[i]]]);
        }
    }
    std::array<char, 5> itoc{'A', 'A', 'C', 'G', 'T'};
    for(const auto& count : counts){
        g0.push_back(itoc[std::distance(count.cbegin(), std::max_element(count.cbegin()+1, count.cend()))]);
    }
    return g0;
}

std::vector<SparseData> to_sparse_data(const std::vector<std::string>& data){
    std::vector<SparseData> sparseData;
    sparseData.reserve(data.size());
    auto lookup = lookupTable();
    auto seq0 = get_reference_expression(data);
    std::size_t count{0};
    for(const auto& seq : data){
        sparseData.emplace_back();
        auto &d = sparseData.back();
        for(std::size_t i=0; i<seq.size(); ++i){
            if(seq0[i] != seq[i]){
                d.push_back(i);
                d.push_back(lookup[seq[i]]);
            }
        }
    }
    return sparseData;
}

std::vector<std::vector<GeneBlock>> to_dense_data(const std::vector<std::string>& data){
    std::vector<std::vector<GeneBlock>> dense;
    dense.reserve(data.size());
    for (const auto &str : data) {
        dense.push_back(from_string(str));
    }
    return dense;
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
