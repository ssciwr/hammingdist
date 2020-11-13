#include"hamming.hh"
#include"hamming_impl.hh"

#include<array>
#include<fstream>
#include<iostream>
#include<sstream>
#include<string>
#include<vector>
#ifdef HAMMING_WITH_OPENMP
#include<omp.h>
#endif

DataSet::DataSet(const std::vector<std::vector<GeneBlock>>& data_)
  : nsamples(data_.size())
  , data(data_)
  , result((nsamples - 1) * nsamples / 2, 0)
{
  validate_data(data);
#ifdef HAMMING_WITH_OPENMP
  #pragma omp parallel for
#endif
  for(std::size_t i=0; i<nsamples; ++i){
    std::size_t offset{i * (i - 1) / 2};
    for(std::size_t j=0; j<i; ++j)
      result[offset + j] = distance(data[i], data[j]);
  }
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
  while(count < n && !stream.eof())
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
