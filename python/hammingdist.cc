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

  py::class_<DataSet<DefaultDistIntType>>(m, "DataSet")
      .def("dump", &DataSet<DefaultDistIntType>::dump,
           "Dump distances matrix in csv format")
      .def("dump_lower_triangular",
           &DataSet<DefaultDistIntType>::dump_lower_triangular,
           "Dump distances matrix in lower triangular format (comma-delimited, "
           "row-major)")
      .def("dump_sequence_indices",
           &DataSet<DefaultDistIntType>::dump_sequence_indices,
           "Dump row index in distances matrix for each input sequence")
      .def("__getitem__", &DataSet<DefaultDistIntType>::operator[])
      .def_readonly("_distances", &DataSet<DefaultDistIntType>::result);

  py::class_<DataSet<uint16_t>>(m, "DataSetLarge")
      .def("dump", &DataSet<uint16_t>::dump,
           "Dump distances matrix in csv format")
      .def("dump_lower_triangular", &DataSet<uint16_t>::dump_lower_triangular,
           "Dump distances matrix in lower triangular format (comma-delimited, "
           "row-major)")
      .def("dump_sequence_indices", &DataSet<uint16_t>::dump_sequence_indices,
           "Dump row index in distances matrix for each input sequence")
      .def("__getitem__", &DataSet<uint16_t>::operator[])
      .def_readonly("_distances", &DataSet<uint16_t>::result);

  m.def("from_stringlist", &from_stringlist,
        "Creates a dataset from a list of strings");
  m.def("from_csv", &from_csv,
        "Creates a dataset by reading already computed distances from csv "
        "(full matrix expected)");
  m.def("from_fasta", &from_fasta<uint8_t>, py::arg("filename"),
        py::arg("include_x") = false, py::arg("remove_duplicates") = false,
        py::arg("n") = 0,
        "Creates a dataset by reading from a fasta file (assuming all "
        "sequences have equal length). Maximum value of an element in the "
        "distances matrix: 255. Distances that would have been larger than "
        "this value saturate at 255 - to support genomes with larger distances "
        "than this see `from_fasta_large` instead.");
  m.def("from_fasta_large", &from_fasta<uint16_t>, py::arg("filename"),
        py::arg("include_x") = false, py::arg("remove_duplicates") = false,
        py::arg("n") = 0,
        "Creates a dataset by reading from a fasta file (assuming all "
        "sequences have equal length). Maximum value of an element in the "
        "distances matrix: 65535");
  m.def("from_lower_triangular", &from_lower_triangular,
        "Creates a dataset by reading already computed distances from lower "
        "triangular format");
  m.def("distance", &distance, py::arg("seq0"), py::arg("seq1"),
        py::arg("include_x") = false,
        "Calculate the distance between seq0 and seq1");
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
  m.def(
      "fasta_sequence_indices",
      [](const std::string &fasta_file, std::size_t n) {
        return as_pyarray(fasta_sequence_indices(fasta_file, n));
      },
      py::arg("fasta_file"), py::arg("n") = 0,
      "Returns the same output as dump_sequence_indices() but without "
      "constructing the distances matrix."
      "For each genome in the input fasta file it gives the index of the "
      "corresponding row in the distances matrix which excludes duplicates");
}

} // namespace hamming
