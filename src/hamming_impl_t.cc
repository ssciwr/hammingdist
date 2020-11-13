#include "hamming_impl.hh"
#include "tests.hh"

using namespace hamming;

TEST_CASE("distance_cpp() returns zero for identical vectors",
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
    REQUIRE(distance_cpp(g1, g1) == 0);
  }
}

TEST_CASE("distance_cpp() returns n for n A's and n G's",
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
    REQUIRE(distance_cpp(g1, g2) == n);
  }
}

TEST_CASE("distance_sparse() returns zero for identical vectors",
          "[impl][distance][sparse]") {
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
    auto g1{make_test_string(n, gen)};
    auto sparse = to_sparse_data({g1, g1});
    REQUIRE(distance_sparse(sparse[0], sparse[1]) == 0);
  }
}

TEST_CASE("distance_sparse() returns n for n A's and n G's",
          "[impl][distance][sparse]") {
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
    auto g1 = std::string(n, 'A');
    auto g2 = std::string(n, 'G');
    auto sparse = to_sparse_data({g1, g2});
    REQUIRE(distance_sparse(sparse[0], sparse[1]) == n);
  }
}

TEST_CASE("distance_sparse() returns same as distance_cpp() for random vectors",
          "[impl][distance][sparse]") {
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
    auto s1{make_test_string(n, gen)};
    auto s2{make_test_string(n, gen)};
    auto g1{from_string(s1)};
    auto g2{from_string(s2)};
    auto sparse = to_sparse_data({s1, s2});
    REQUIRE(distance_sparse(sparse[0], sparse[1]) == distance_cpp(g1, g2));
  }
}

