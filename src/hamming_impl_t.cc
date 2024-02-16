#include "hamming/hamming_impl.hh"
#include "tests.hh"

using namespace hamming;

TEST_CASE("distance_cpp() returns zero for identical vectors of valid chars",
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
    for (int max_dist : {0, 1, 2, 11, 999, 9876544}) {
      CAPTURE(n);
      CAPTURE(max_dist);
      auto g1{make_gene_vector(n, gen)};
      REQUIRE(distance_cpp(g1, g1, max_dist) == 0);
    }
  }
}

TEST_CASE("distance_cpp() returns n for n A's and n G's", "[impl][distance]") {
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
      REQUIRE(distance_cpp(g1, g2, max_dist) == std::min(max_dist, n));
    }
  }
}

TEST_CASE("distance_sparse() returns zero for identical vectors of valid chars",
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
    for (int max_dist : {0, 1, 2, 11, 999, 9876544}) {
      CAPTURE(n);
      CAPTURE(max_dist);
      for (bool include_x : {false, true}) {
        auto g1{make_test_string(n, gen, include_x)};
        CAPTURE(include_x);
        auto sparse = to_sparse_data({g1, g1}, include_x);
        REQUIRE(distance_sparse(sparse[0], sparse[1], max_dist) == 0);
      }
    }
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
    for (int max_dist : {0, 1, 2, 11, 999, 9876544}) {
      CAPTURE(n);
      CAPTURE(max_dist);
      for (char c1 : {'A', 'X'}) {
        for (char c2 : {'G', 'T'}) {
          auto g1 = std::string(n, c1);
          auto g2 = std::string(n, c2);
          for (bool include_x : {false, true}) {
            CAPTURE(include_x);
            auto sparse = to_sparse_data({g1, g2}, include_x);
            REQUIRE(distance_sparse(sparse[0], sparse[1], max_dist) ==
                    std::min(max_dist, n));
          }
        }
      }
    }
  }
}

TEST_CASE("distance_sparse() returns same as distance_cpp() for random vectors "
          "(with no X chars)",
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
    for (int max_dist : {0, 1, 2, 11, 999, 9876544}) {
      CAPTURE(n);
      CAPTURE(max_dist);
      auto s1{make_test_string(n, gen)};
      auto s2{make_test_string(n, gen)};
      auto g1{from_string(s1)};
      auto g2{from_string(s2)};
      for (bool include_x : {false, true}) {
        CAPTURE(include_x);
        auto sparse = to_sparse_data({s1, s2}, include_x);
        REQUIRE(distance_sparse(sparse[0], sparse[1], max_dist) ==
                distance_cpp(g1, g2, max_dist));
      }
    }
  }
}

TEST_CASE("distance_sparse() returns same as distance_cpp() for equal-distance "
          "random vectors with X's",
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
    for (int max_dist : {0, 1, 2, 11, 999, 9876544}) {
      CAPTURE(n);
      CAPTURE(max_dist);
      // calculate distance of strings with no X's
      auto s1{make_test_string(n, gen)};
      auto s2{make_test_string(n, gen)};
      auto g1{from_string(s1)};
      auto g2{from_string(s2)};
      auto dist = distance_cpp(g1, g2, max_dist);
      // replace all A's with X's : sparse distance with X as valid char should
      // be the same
      for (char replace_char : {'A', 'C', 'G', 'T'}) {
        auto s1_x = s1;
        auto s2_x = s2;
        for (auto &c : s1_x) {
          if (c == replace_char) {
            c = 'X';
          }
        }
        for (auto &c : s2_x) {
          if (c == replace_char) {
            c = 'X';
          }
        }
        auto sparse = to_sparse_data({s1_x, s2_x}, true);
        REQUIRE(distance_sparse(sparse[0], sparse[1], max_dist) == dist);
      }
    }
  }
}
