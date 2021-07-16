#ifndef HAMMING_TESTS_HH
#define HAMMING_TESTS_HH

#include <catch2/catch.hpp>
#include <random>
#include <vector>
#include "hamming/hamming.hh"
#include "hamming_impl.hh"

namespace hamming {

std::string make_test_string(int n, std::mt19937 &gen, bool include_x=false);

std::vector<GeneBlock> make_gene_vector(int n, std::mt19937 &gen, bool include_x=false);

}

#endif
