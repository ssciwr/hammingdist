#pragma once

#include "hamming/hamming.hh"
#include "hamming/hamming_impl.hh"
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <random>
#include <vector>

namespace hamming {

std::string make_test_string(int n, std::mt19937 &gen, bool include_x = false);

std::vector<GeneBlock> make_gene_vector(int n, std::mt19937 &gen,
                                        bool include_x = false);

void write_test_fasta(const std::string &filename, int n, std::size_t n_seq,
                      std::mt19937 &gen, bool include_x = false);

} // namespace hamming
