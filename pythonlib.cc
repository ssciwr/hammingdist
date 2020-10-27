#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include<array>
#include<cstdint>
#include<fstream>
#include<iostream>
#include<sstream>
#include<vector>

using Gene = std::uint_fast8_t;

static int distance(const std::vector<Gene>& a, const std::vector<Gene>& b){
  int r{0};
  for (std::size_t i=0; i<a.size(); ++i)
    if(a[i] ^ b[i])
      ++r;
  return r;
}

struct DataSet 
{
  DataSet(const std::vector<std::vector<Gene>>& data_)
    : nsamples(data_.size())
    , data(data_)
    , result((nsamples - 1) * nsamples / 2)
  {
    calculate();
  }

  DataSet(const std::string& filename)
  {
    std::cout << "Not yet implemented" << std::endl;
  }

  void calculate()
  {
    std::size_t pos = 0;
    for(std::size_t i=0; i<nsamples; ++i)
      for(std::size_t j=0; j<i; ++j)
        result[pos++] = distance(data[i], data[j]);
  }

  void dump(const std::string& filename)
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

  int operator[](const std::array<std::size_t, 2>& index) const
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

  std::size_t nsamples;
  std::vector<std::vector<Gene>> data;
  std::vector<int> result;
};

DataSet from_stringlist(const std::vector<std::string>& data)
{
  std::vector<std::vector<Gene>> result(data.size());
  char buffer;

  std::array<Gene, 256> lookup;
  std::fill(lookup.begin(), lookup.end(), 0);
  lookup[std::size_t('A')] = 1 << 1;
  lookup[std::size_t('C')] = 1 << 2;
  lookup[std::size_t('G')] = 1 << 3;
  lookup[std::size_t('T')] = 1 << 4;

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

namespace py = pybind11;

PYBIND11_MODULE(hammingdist,m)
{
  m.doc() = "Small to calculate Hamming distances between gene sequences";

  py::class_<DataSet>(m, "DataSet")
      .def("dump", &DataSet::dump)
      .def("__getitem__", &DataSet::operator[])
      .def_readonly("data", &DataSet::data)
      .def_readonly("distances", &DataSet::result);

  m.def("from_stringlist", &from_stringlist, "Creates a dataset from a list of strings");
  m.def("from_csv", &from_csv, "Creates a dataset by reading already computed distances from csv (full matrix expected)");
}