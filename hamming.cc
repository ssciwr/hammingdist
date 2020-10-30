#include"hamming.hh"

#include<array>
#include<fstream>
#include<iostream>
#include<sstream>
#include<string>
#include<vector>
#include<limits>
#include<stdexcept>
#ifdef HAMMING_WITH_OPENMP
#include<omp.h>
#endif
// 4-bit representation of gene:
constexpr std::size_t n_bits_per_gene{4};
constexpr GeneBlock mask_gene0{0x0f};
constexpr GeneBlock mask_gene1{0xf0};

// bit meaning:
// 1111: '-'
// 0001: 'A'
// 0010: 'C'
// 0100: 'G'
// 1000: 'T'
// 0000: invalid
static std::array<GeneBlock, 256> lookupTable()
{
  std::array<GeneBlock, 256> lookup{};
  lookup[std::size_t('-')] = 0xff;
  lookup[std::size_t('A')] = 1 | (1 << n_bits_per_gene);
  lookup[std::size_t('C')] = (1 << 1) | (1 << (n_bits_per_gene+1));
  lookup[std::size_t('G')] = (1 << 2) | (1 << (n_bits_per_gene+2));
  lookup[std::size_t('T')] = (1 << 3) | (1 << (n_bits_per_gene+3));

  return lookup;
}

static int distance(const std::vector<GeneBlock>& a, const std::vector<GeneBlock>& b){
  int r{0};
  for (std::size_t i=0; i<a.size(); ++i) {
    auto c{static_cast<GeneBlock>(a[i] & b[i])};
    r += static_cast<int>((c & mask_gene0) == 0);
    r += static_cast<int>((c & mask_gene1) == 0);
  }
  return r;
}

static void validate_data(const std::vector<std::vector<GeneBlock>>& data){
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

DataSet::DataSet(const std::vector<std::vector<GeneBlock>>& data_)
  : nsamples(data_.size())
  , data(data_)
  , result((nsamples - 1) * nsamples / 2, 0)
{
  validate_data(data);
  std::size_t pos = 0;
#ifdef HAMMING_WITH_OPENMP
  #pragma omp parallel for
#endif
  for(std::size_t i=0; i<nsamples; ++i)
    for(std::size_t j=0; j<i; ++j)
      result[pos++] = distance(data[i], data[j]);
}

DataSet::DataSet(const std::string& filename)
{
  std::cout << "Not yet implemented" << std::endl;
}

void DataSet::dump(const std::string& filename)
{
  std::ofstream stream(filename);
  for(std::size_t i=0; i<nsamples; ++i)
  {
    for(std::size_t j=0; j<nsamples; ++j)
    {
      stream << (*this)[{i, j}];
      if (j != nsamples - 1)
        stream << ", ";
    }
    stream << std::endl;
  }
}

int DataSet::operator[](const std::array<std::size_t, 2>& index) const
{
  auto i = index[0];
  auto j = index[1];
  if (i < j)
    return result[j * (j - 1) / 2 + i];
  if (i > j)
    return result[i * (i - 1) / 2 + j];
  // This is a diagonal entry
  return 0;
}

static std::vector<GeneBlock> from_string(const std::string& str){
  std::vector<GeneBlock> r;
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

DataSet from_stringlist(const std::vector<std::string> &data) {
  std::vector<std::vector<GeneBlock>> result;
  result.reserve(data.size());
  for (const auto &str : data) {
      result.push_back(from_string(str));
  }
  return DataSet(result);
}

DataSet from_csv(const std::string& filename)
{
  return DataSet(filename);
}

DataSet from_fasta(const std::string& filename, std::size_t n)
{
  std::vector<std::vector<GeneBlock>> m;
  m.reserve(n);
  // Initializing the stream
  std::ifstream stream(filename);
  std::size_t count = 0;
  std::string line;
  // skip first header
  std::getline(stream, line);
  while(count < n)
  {
    std::string seq;
    while(std::getline(stream, line) && line[0] != '>')
    {
      seq.append(line);
    }
    m.push_back(from_string(seq));
    ++count;
  }
  return DataSet(m);
}
