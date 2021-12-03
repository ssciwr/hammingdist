#ifndef _HAMMING_HH
#define _HAMMING_HH

#include "hamming/hamming_types.hh"
#include <string>
#include <vector>

namespace hamming {

DataSet from_stringlist(std::vector<std::string> &);
DataSet from_csv(const std::string &);
DataSet from_fasta(const std::string &, bool include_x = false,
                   bool remove_duplicates = false, std::size_t n = 0);
DataSet from_lower_triangular(const std::string &);

ReferenceDistIntType distance(const std::string &seq0, const std::string &seq1,
                              bool include_x = false);
std::vector<ReferenceDistIntType>
fasta_reference_distances(const std::string &reference_sequence,
                          const std::string &fasta_file,
                          bool include_x = false);

} // namespace hamming

#endif
