#define CATCH_CONFIG_MAIN

#include "tests.hh"

namespace hamming {

std::string make_test_string(int n, std::mt19937 &gen) {
  std::string s;
  std::uniform_int_distribution<> distrib(0, 4);
  std::array<char, 5> c{'A', 'C', 'G', 'T', '-'};
  s.reserve(n);
  for (int i = 0; i < n; ++i) {
    s.push_back(c[distrib(gen)]);
  }
  return s;
}

std::vector<GeneBlock> make_gene_vector(int n, std::mt19937 &gen) {
  return from_string(make_test_string(n, gen));
}

}
