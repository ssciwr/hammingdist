#include <memory>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "hamming/hamming.hh"

namespace py = pybind11;

// helper function to avoid making a copy when returning a py::array_t
// author: https://github.com/YannickJadoul
// source: https://github.com/pybind/pybind11/issues/1042#issuecomment-642215028
template <typename Sequence>
inline py::array_t<typename Sequence::value_type> as_pyarray(Sequence &&seq) {
  auto size = seq.size();
  auto data = seq.data();
  std::unique_ptr<Sequence> seq_ptr =
      std::make_unique<Sequence>(std::move(seq));
  auto capsule = py::capsule(seq_ptr.get(), [](void *p) {
    std::unique_ptr<Sequence>(reinterpret_cast<Sequence *>(p));
  });
  seq_ptr.release();
  return py::array(size, data, capsule);
}

namespace hamming {

PYBIND11_MODULE(hammingdist, m) {
  m.doc() = "Small tool to calculate Hamming distances between gene sequences";

  py::class_<DataSet>(m, "DataSet")
      .def("dump", &DataSet::dump, "Dump distances matrix in csv format")
      .def("dump_lower_triangular", &DataSet::dump_lower_triangular,
           "Dump distances matrix in lower triangular format (comma-delimited, "
           "row-major)")
      .def("dump_sequence_indices", &DataSet::dump_sequence_indices,
           "Dump row index in distances matrix for each input sequence")
      .def("__getitem__", &DataSet::operator[])
      .def_readonly("_distances", &DataSet::result);

  m.def("from_stringlist", &from_stringlist,
        "Creates a dataset from a list of strings");
  m.def("from_csv", &from_csv,
        "Creates a dataset by reading already computed distances from csv "
        "(full matrix expected)");
  m.def("from_fasta", &from_fasta, py::arg("filename"),
        py::arg("include_x") = false, py::arg("remove_duplicates") = false,
        py::arg("n") = 0,
        "Creates a dataset by reading from a fasta file (assuming all "
        "sequences have equal length)");
  m.def("from_lower_triangular", &from_lower_triangular,
        "Creates a dataset by reading already computed distances from lower "
        "triangular format");
  m.def(
      "fasta_reference_distances",
      [](const std::string &reference_sequence, const std::string &fasta_file,
         bool include_x) {
        return as_pyarray(fasta_reference_distances(reference_sequence,
                                                    fasta_file, include_x));
      },
      py::arg("reference_sequence"), py::arg("fasta_file"),
      py::arg("include_x") = false,
      "Calculates the distance of each sequence in the fasta file from the "
      "supplied reference sequence");
}

} // namespace hamming
