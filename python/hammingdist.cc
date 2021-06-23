#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "hamming/hamming.hh"

namespace py = pybind11;

namespace hamming {

PYBIND11_MODULE(hammingdist,m)
{
  m.doc() = "Small tool to calculate Hamming distances between gene sequences";

  py::class_<DataSet>(m, "DataSet")
      .def("dump", &DataSet::dump, "Dump distances matrix in csv format")
      .def("dump_lower_triangular", &DataSet::dump_lower_triangular, "Dump distances matrix in lower triangular format (comma-delimited, row-major)")
      .def("__getitem__", &DataSet::operator[])
      .def_readonly("_distances", &DataSet::result);

  m.def("from_stringlist", &from_stringlist, "Creates a dataset from a list of strings");
  m.def("from_csv", &from_csv, "Creates a dataset by reading already computed distances from csv (full matrix expected)");
  m.def("from_fasta", &from_fasta, py::arg("filename"), py::arg("n") = 0, "Creates a dataset by reading from a fasta file (assuming all sequences have equal length)");
  m.def("from_lower_triangular", &from_lower_triangular, "Creates a dataset by reading already computed distances from lower triangular format");
}

}
