#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "hamming.hh"

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
  m.def("from_fasta", &from_fasta, "Creates a dataset by reading from a fasta file (assuming all sequences have equal length)");
}