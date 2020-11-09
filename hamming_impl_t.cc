#include "hamming_impl.hh"
#include <catch2/catch.hpp>
#include <random>
#include <string>
#include <vector>

static std::vector<GeneBlock> make_gene_vector(int n, std::mt19937 &gen) {
  std::string s;
  std::uniform_int_distribution<> distrib(0, 4);
  std::array<char, 5> c{'A', 'C', 'G', 'T', '-'};
  s.reserve(n);
  for (int i = 0; i < n; ++i) {
    s.push_back(c[distrib(gen)]);
  }
  return from_string(s);
}

TEST_CASE("distance() and distance_simd() give zero for identical vectors",
          "[impl][distance]") {
  std::mt19937 gen(12345);
  for (int n :
       {1,       2,      3,      4,      5,      6,      7,      8,
        9,       10,     11,     12,     13,     14,     15,     16,
        17,      18,     19,     20,     31,     32,     33,     63,
        64,      65,     127,    128,    129,    254,    255,    256,
        256,     511,    512,    513,    1023,   1024,   1025,   2047,
        2048,    2049,   4095,   4096,   4097,   8191,   8192,   8193,
        32767,   32768,  32769,  65535,  65536,  65537,  131071, 131072,
        131073,  262143, 262144, 262145, 524287, 524288, 524289, 1048575,
        1048576, 1048577}) {
    CAPTURE(n);
    auto g1{make_gene_vector(n, gen)};
    REQUIRE(distance(g1, g1) == 0);
    REQUIRE(distance_sse2(g1, g1) == 0);
  }
}

TEST_CASE("distance() and distance_simd() give n for n A's and n G's",
          "[impl][distance]") {
  for (int n :
       {1,       2,      3,      4,      5,      6,      7,      8,
        9,       10,     11,     12,     13,     14,     15,     16,
        17,      18,     19,     20,     31,     32,     33,     63,
        64,      65,     127,    128,    129,    254,    255,    256,
        256,     511,    512,    513,    1023,   1024,   1025,   2047,
        2048,    2049,   4095,   4096,   4097,   8191,   8192,   8193,
        32767,   32768,  32769,  65535,  65536,  65537,  131071, 131072,
        131073,  262143, 262144, 262145, 524287, 524288, 524289, 1048575,
        1048576, 1048577}) {
    CAPTURE(n);
    auto g1 = from_string(std::string(n, 'A'));
    auto g2 = from_string(std::string(n, 'G'));
    REQUIRE(distance(g1, g2) == n);
    REQUIRE(distance_sse2(g1, g2) == n);
  }
}

TEST_CASE("distance() and distance_simd() give the same results for random vectors",
          "[impl][distance]") {
  std::mt19937 gen(12345);
  for (int n :
       {1,       2,      3,      4,      5,      6,      7,      8,
        9,       10,     11,     12,     13,     14,     15,     16,
        17,      18,     19,     20,     31,     32,     33,     63,
        64,      65,     127,    128,    129,    254,    255,    256,
        256,     511,    512,    513,    1023,   1024,   1025,   2047,
        2048,    2049,   4095,   4096,   4097,   8191,   8192,   8193,
        32767,   32768,  32769,  65535,  65536,  65537,  131071, 131072,
        131073,  262143, 262144, 262145, 524287, 524288, 524289, 1048575,
        1048576, 1048577}) {
    CAPTURE(n);
    auto g1{make_gene_vector(n, gen)};
    auto g2{make_gene_vector(n, gen)};
    REQUIRE(distance(g1, g2) == distance_sse2(g1, g2));
  }
}
