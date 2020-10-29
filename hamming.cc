#include"hamming.hh"

#include<array>
#include<fstream>
#include<iostream>
#include<sstream>
#include<string>
#include<vector>
#include<limits>

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

DataSet::DataSet(const std::vector<std::vector<GeneBlock>>& data_)
  : nsamples(data_.size())
  , data(data_)
  , result((nsamples - 1) * nsamples / 2, 0)
{
  std::size_t pos = 0;
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

DataSet from_stringlist(const std::vector<std::string> &data) {
  std::vector<std::vector<GeneBlock>> result;
  result.reserve(data.size());
  auto lookup = lookupTable();
  for (const auto &str : data) {
    std::size_t n_full_blocks{str.size()/2};
    result.emplace_back();
    auto &r = result.back();
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
  }
  return DataSet(result);
}

DataSet from_csv(const std::string& filename)
{
  return DataSet(filename);
}

DataSet from_fasta(const std::string& filename, std::size_t n)
{
  // Determine the length of the sequence in the given fasta file
  unsigned int l = 0;
  std::ifstream istream(filename);
  std::string line;
  // Consume the header line
  std::getline(istream, line);
  while(true)
  {
    std::getline(istream, line);
    if(line.substr(0, 1) == ">")
      break;
    else
      l += line.size();
  }
  istream.close();

  std::vector<std::vector<GeneBlock>> m;
  m.reserve(n);

  // Initializing the stream
  std::ifstream stream(filename);
  std::size_t count = 0;
  unsigned char buffer;
  auto lookup = lookupTable();

  // The while-condition consumes the header line
  while((count < n) && (std::getline(stream, line)))
  {
    m.emplace_back();
    auto &v = m.back();
    v.resize(l);
    std::size_t pos = 0;

    while((pos < l) && (std::getline(stream, line)))
    {
      std::istringstream iss(line);
      while(iss >> buffer)
        v[pos++] = lookup[std::size_t(buffer)];
    }

    ++count;
  }

  return DataSet(m);
}
