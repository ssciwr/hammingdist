#include"hamming/hamming.hh"
#include"hamming_impl.hh"

#include<array>
#include<fstream>
#include<iostream>
#include<sstream>
#include<string>
#include<vector>

namespace hamming {

DataSet::DataSet(const std::vector<std::string>& data_)
  : nsamples(data_.size())
{
  validate_data(data_);
  result = distances(data_);
}

DataSet::DataSet(const std::string& filename)
{
  std::ifstream stream(filename);
  std::string line;
  nsamples = 0;

  while(std::getline(stream, line))
  {
    std::istringstream s(line);
    std::string d;
    for(std::size_t j=0; j<nsamples; ++j)
    {
      std::getline(s, d, ',');
      result.push_back(std::stoi(d));
    }
    ++nsamples;
  }
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
  return DataSet(data);
}

DataSet from_csv(const std::string& filename)
{
  return DataSet(filename);
}

DataSet from_fasta(const std::string& filename, std::size_t n)
{
  std::vector<std::string> data;
  data.reserve(n);
  // Initializing the stream
  std::ifstream stream(filename);
  std::size_t count = 0;
  std::string line;
  // skip first header
  std::getline(stream, line);
  while(count < n && !stream.eof())
  {
    data.emplace_back();
    auto& seq = data.back();
    while(std::getline(stream, line) && line[0] != '>')
    {
      seq.append(line);
    }
    ++count;
  }
  return DataSet(data);
}

}
