#include"hamming/hamming.hh"
#include"hamming_impl.hh"

#include<array>
#include<algorithm>
#include<cmath>
#include<fstream>
#include<iostream>
#include<limits>
#include<sstream>
#include<string>
#include<vector>

namespace hamming {

DataSet::DataSet(std::vector<std::string>& data_, bool clear_input_data)
  : nsamples(data_.size())
{
  validate_data(data_);
  result = distances(data_, clear_input_data);
}

DataSet::DataSet(const std::string& filename)
{
  // Determine correct dataset size
  std::ifstream stream(filename);
  std::string line;
  nsamples = std::count(std::istreambuf_iterator<char>(stream),
                        std::istreambuf_iterator<char>(), '\n');
  result.resize(nsamples * (nsamples + 1) / 2);

  // Read the data
  stream = std::ifstream(filename);
  std::size_t i{0};
  std::size_t current{0};
  while(std::getline(stream, line))
  {
    std::istringstream s(line);
    std::string d;
    for(std::size_t j=0; j<current; ++j)
    {
      std::getline(s, d, ',');
      result[i++] = std::stoi(d);
    }
    ++current;
  }
}

static std::size_t uint_sqrt(std::size_t x){
  return static_cast<std::size_t>(std::round(std::sqrt(x)));
}

DataSet::DataSet(std::vector<DistIntType>&& distances) : result{std::move(distances)}
{
  // infer n from number of lower triangular matrix elements = n(n-1)/2
  nsamples = (uint_sqrt(8*result.size()+1)+1)/2;
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

void DataSet::dump_lower_triangular(const std::string& filename)
{
  std::ofstream stream(filename);
  std::size_t k = 0;
  for(std::size_t i=1; i<nsamples; ++i)
  {
    for(std::size_t j=0; j+1<i; ++j)
    {
      stream << static_cast<int>(result[k++]) << ",";
    }
    stream << static_cast<int>(result[k++]) << "\n";
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

DataSet from_stringlist(std::vector<std::string> &data) {
  return DataSet(data);
}

DataSet from_csv(const std::string& filename)
{
  return DataSet(filename);
}

DataSet from_lower_triangular(const std::string& filename)
{
  std::vector<DistIntType> distances;
  std::ifstream stream(filename);
  std::string line;
  while(std::getline(stream, line))
  {
    std::istringstream s(line);
    std::string d;
    while(s.good())
    {
      std::getline(s, d, ',');
      distances.push_back(safe_int_cast(std::stoi(d)));
    }
  }
  return DataSet(std::move(distances));
}

DataSet from_fasta(const std::string& filename, std::size_t n)
{
  std::vector<std::string> data;
  data.reserve(n);
  if (n == 0) {
    n = std::numeric_limits<std::size_t>::max();
    data.reserve(65536);
  }
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
  return DataSet(data, true);
}

}
