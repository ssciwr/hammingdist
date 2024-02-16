#include "hamming/distance_neon.hh"
#include "tests.hh"

using namespace hamming;

TEST_CASE("distance_neon() returns all return zero for identical vectors",
          "[impl][distance][neon]") {
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
    for (int max_dist : {0, 1, 2, 11, 999, 9876544}) {
      CAPTURE(n);
      CAPTURE(max_dist);
      auto g1{make_gene_vector(n, gen)};
      REQUIRE(distance_neon(g1, g1, max_dist) == 0);
    }
  }
}

TEST_CASE("distance_neon() all return n for n A's and n G's",
          "[impl][distance][neon]") {
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
    for (int max_dist : {0, 1, 2, 11, 999, 9876544}) {
      CAPTURE(n);
      CAPTURE(max_dist);
      auto g1 = from_string(std::string(n, 'A'));
      auto g2 = from_string(std::string(n, 'G'));
      REQUIRE(distance_neon(g1, g2, max_dist) == std::min(max_dist, n));
    }
  }
}

TEST_CASE("distance_neon() returns same as distance_cpp() for random vectors",
          "[impl][distance][neon]") {
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
    for (int max_dist : {0, 1, 2, 11, 999, 9876544}) {
      CAPTURE(n);
      CAPTURE(max_dist);
      auto g1{make_gene_vector(n, gen)};
      auto g2{make_gene_vector(n, gen)};
      REQUIRE(distance_neon(g1, g2, max_dist) ==
              distance_cpp(g1, g2, max_dist));
    }
  }
}
