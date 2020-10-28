#include"hamming.hh"

#include<array>
#include<fstream>
#include<iostream>
#include<sstream>
#include<string>
#include<vector>

static std::array<Gene, 256> lookupTable()
{
  std::array<Gene, 256> lookup;
  std::fill(lookup.begin(), lookup.end(), 0xFF);
  lookup[std::size_t('A')] = 1 << 1;
  lookup[std::size_t('C')] = 1 << 2;
  lookup[std::size_t('G')] = 1 << 3;
  lookup[std::size_t('T')] = 1 << 4;

  return lookup;
}

static int distance(const std::vector<Gene>& a, const std::vector<Gene>& b){
  int r{0};
  for (std::size_t i=0; i<a.size(); ++i)
      r += static_cast<int>((a[i] & b[i]) == 0);
  return r;
}

DataSet::DataSet(const std::vector<std::vector<Gene>>& data_)
  : nsamples(data_.size())
  , data(data_)
  , result((nsamples - 1) * nsamples / 2)
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

DataSet from_stringlist(const std::vector<std::string>& data)
{
  std::vector<std::vector<Gene>> result(data.size());
  char buffer;

  auto lookup = lookupTable();
  for(std::size_t i =0; i < data.size(); ++i)
  {
    result[i].resize(data[i].size());
    std::size_t pos = 0;
    std::istringstream iss(data[i]);
    while(iss >> buffer)
      result[i][pos++] = lookup[std::size_t(buffer)];
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

  std::vector<std::vector<Gene>> m;
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
