#include "hamming/hamming.hh"
#include "tests.hh"
#include <cstdio>
#include <fstream>
#include <string>

using namespace hamming;

static constexpr std::array<char, 4> valid_chars{'A', 'C', 'G', 'T'};
static constexpr std::array<char, 6> invalid_chars{' ', 'N', '*',
                                                   '?', 'a', '.'};
static std::vector<bool> valid_use_gpu_values() {
  if (cuda_gpu_available()) {
    return {false, true};
  }
  return {false};
}

static int dist1(char c1, char c2) {
  std::vector<std::string> v{std::string{c1}, std::string{c2}};
  return from_stringlist(v)[{0, 1}];
}

static int dist2(char c1, char c2) {
  return static_cast<int>(distance(std::string{c1}, std::string{c2}));
}

TEST_CASE("distance between two equal valid characters is 0", "[distance]") {
  for (auto c : valid_chars) {
    CAPTURE(c);
    REQUIRE(dist1(c, c) == 0);
    REQUIRE(dist2(c, c) == 0);
  }
}

TEST_CASE("distance between any valid char & '-' is 0", "[distance]") {
  for (auto c : valid_chars) {
    CAPTURE(c);
    REQUIRE(dist1(c, '-') == 0);
    REQUIRE(dist1('-', c) == 0);
    REQUIRE(dist2(c, '-') == 0);
    REQUIRE(dist2('-', c) == 0);
  }
}

TEST_CASE("distance between two different valid characters is 1",
          "[distance]") {
  for (auto c1 : valid_chars) {
    for (auto c2 : valid_chars) {
      if (c1 != c2) {
        CAPTURE(c1);
        CAPTURE(c2);
        REQUIRE(dist1(c1, c2) == 1);
        REQUIRE(dist2(c1, c2) == 1);
      }
    }
  }
}

TEST_CASE("distance between valid & invalid characters is 1", "[distance]") {
  for (auto c1 : valid_chars) {
    for (auto c2 : invalid_chars) {
      CAPTURE(c1);
      CAPTURE(c2);
      CAPTURE((int)c2);
      REQUIRE(dist1(c1, c2) == 1);
      REQUIRE(dist2(c1, c2) == 1);
    }
  }
}

TEST_CASE("distance between two invalid characters is 1", "[distance]") {
  for (auto c1 : invalid_chars) {
    for (auto c2 : invalid_chars) {
      CAPTURE(c1);
      CAPTURE(c2);
      REQUIRE(dist1(c1, c2) == 1);
      REQUIRE(dist2(c1, c2) == 1);
    }
  }
}

TEST_CASE("distance between any invalid char & '-' is 1", "[distance]") {
  for (auto c : invalid_chars) {
    CAPTURE(c);
    REQUIRE(dist1(c, '-') == 1);
    REQUIRE(dist1('-', c) == 1);
    REQUIRE(dist2(c, '-') == 1);
    REQUIRE(dist2('-', c) == 1);
  }
}

TEST_CASE("two expressions with distance 2", "[hamming]") {
  std::vector<std::vector<std::string>> expr;
  expr.push_back({{"AC"}, {"CA"}});
  expr.push_back({{"AC"}, {"TG"}});
  expr.push_back({{"ACG"}, {"AGT"}});
  expr.push_back({{"ACG"}, {"-TT"}});
  expr.push_back({{"ACG"}, {"T-T"}});
  expr.push_back({{"ACG"}, {"TA-"}});
  expr.push_back({{"ACG"}, {"CCC"}});
  expr.push_back({{"ACGT"}, {"AGGG"}});
  expr.push_back({{"ACGTGTCGTGTCGACGTGTCG"}, {"ACGTGTCGTTTCGACGAGTCG"}});
  expr.push_back({{"ACGTGTCGTGTCGACGTGTCGT"}, {"ACGTGTCGTTTCGACGAGTCGT"}});
  expr.push_back({{"ACGTGTCGTGTCGACGTGTCGT-"}, {"ACGTGTCGTTTCGACGAGTCGTA"}});
  expr.push_back({{"ACGTGTCGTGTCGACGTGTCG----------"},
                  {"ACGTGTCGTTTCGACGAGTCGGGG-------"}});
  expr.push_back({{"ACGTGTCGTGTCGACGTGTCG----------A"},
                  {"ACGTGTCGTTTCGACGAGTCGGGG-------A"}});
  expr.push_back({{"ACGTGTCGTGTCGACGTGTCG----------AG"},
                  {"ACGTGTCGTTTCGACGAGTCGGGG-------AG"}});
  expr.push_back({{"ACGTGTCGTGTCGACGTGTCG---ACGTGTCGTGTCGACGTGTCG---"
                   "ACGTGTCGTGTCGACGTGTCG---"},
                  {"ACGTGTCGTTTCGACGAGTCGGGGACGTGTCGTGTCGACGTGTCG---"
                   "ACGTGTCGTGTCGACGTGTCG---"}});
  for (auto &v : expr) {
    auto d = from_stringlist(v);
    REQUIRE(d[{0, 0}] == 0);
    REQUIRE(distance(v[0], v[0]) == 0);
    REQUIRE(d[{0, 1}] == 2);
    REQUIRE(distance(v[0], v[1]) == 2);
    REQUIRE(d[{1, 0}] == 2);
    REQUIRE(distance(v[1], v[0]) == 2);
    REQUIRE(d[{1, 1}] == 0);
    REQUIRE(distance(v[1], v[1]) == 0);
  }
}

TEST_CASE("from_fasta single line sequences", "[hamming]") {
  char tmp_file_name[L_tmpnam];
  REQUIRE(std::tmpnam(tmp_file_name) != nullptr);
  CAPTURE(tmp_file_name);
  std::ofstream of(tmp_file_name);
  of << ">seq0\n";
  of << "ACGTGTCGTGTCGACGTGTCG\n";
  of << ">seq1\n";
  of << "ACGTGTCGTTTCGACGAGTCG\n";
  of.close();
  for (auto use_gpu : valid_use_gpu_values()) {
    for (bool remove_duplicates : {false, true}) {
      for (bool include_x : {false, true}) {
        CAPTURE(include_x);
        CAPTURE(remove_duplicates);
        CAPTURE(use_gpu);
        if (use_gpu && include_x) {
          // include_x cannot be used with use_gpu
          REQUIRE_THROWS(from_fasta<uint8_t>(tmp_file_name, include_x,
                                             remove_duplicates, 1, use_gpu));
        } else {
          for (int n : {0, 2, 3, 8}) {
            auto d = from_fasta<uint8_t>(tmp_file_name, include_x,
                                         remove_duplicates, n, use_gpu);
            REQUIRE(d[{0, 0}] == 0);
            REQUIRE(d[{0, 1}] == 2);
            REQUIRE(d[{1, 0}] == 2);
            REQUIRE(d[{1, 1}] == 0);
          }
        }
      }
    }
  }
  std::remove(tmp_file_name);
}

TEST_CASE("from_fasta single line sequences with duplicates", "[hamming]") {
  char tmp_file_name[L_tmpnam];
  REQUIRE(std::tmpnam(tmp_file_name) != nullptr);
  CAPTURE(tmp_file_name);
  std::ofstream of(tmp_file_name);
  of << ">seq0\n";
  of << "ACGTGTCGTGTCGACGTGTCG\n";
  of << ">seq1\n";
  of << "ACGTGTCGTTTCGACGAGTCG\n";
  of << ">seq2\n";
  of << "ACGTGTCGTTTCGACGAGTCG\n";
  of << ">seq3\n";
  of << "ACGTGTCGTTTCGACGAGTCG\n";
  of << ">seq4\n";
  of << "ACGTGTCGTGTCGACGTGTCG\n";
  of << ">seq5\n";
  of << "ACGTGTCGTATCGACGTGTCG\n";
  of.close();
  std::vector<std::size_t> sequence_indices{0, 1, 1, 1, 0, 2};
  for (auto use_gpu : valid_use_gpu_values()) {
    for (int n : {0, 6, 22}) {
      for (bool include_x : {false, true}) {
        CAPTURE(include_x);
        CAPTURE(use_gpu);
        if (use_gpu && include_x) {
          // include_x cannot be used with use_gpu
          REQUIRE_THROWS(
              from_fasta<uint8_t>(tmp_file_name, include_x, true, n, use_gpu));
        } else {
          auto d =
              from_fasta<uint8_t>(tmp_file_name, include_x, true, n, use_gpu);
          REQUIRE(d.nsamples == 3);
          REQUIRE(d.sequence_indices == sequence_indices);
          REQUIRE(d[{0, 0}] == 0);
          REQUIRE(d[{0, 1}] == 2);
          REQUIRE(d[{0, 2}] == 1);
          REQUIRE(d[{1, 0}] == 2);
          REQUIRE(d[{1, 1}] == 0);
          REQUIRE(d[{1, 2}] == 2);
          REQUIRE(d[{2, 0}] == 1);
          REQUIRE(d[{2, 1}] == 2);
          REQUIRE(d[{2, 2}] == 0);
        }
      }
    }
  }
  std::remove(tmp_file_name);
}

TEST_CASE("from_fasta multi-line sequences", "[hamming]") {
  char tmp_file_name[L_tmpnam];
  REQUIRE(std::tmpnam(tmp_file_name) != nullptr);
  CAPTURE(tmp_file_name);
  std::ofstream of(tmp_file_name);
  of << ">seq0\n";
  of << "ACGTGTCGTGTCGACGTGTCG\n";
  of << "ACGTGTCGTGTCGACGTGTCG\n";
  of << "ACGTGTCGTGTCG\n";
  of << ">seq1\n";
  of << "ACGTGTCGTTTCGACGAGTCG\n";
  of << "ACGTGTCGTGTCGACGTGTCG\n";
  of << "ACGTGTCGTGTCG\n";
  of.close();
  for (auto use_gpu : valid_use_gpu_values()) {
    for (bool remove_duplicates : {false, true}) {
      for (bool include_x : {false, true}) {
        CAPTURE(include_x);
        CAPTURE(use_gpu);
        if (use_gpu && include_x) {
          // include_x cannot be used with use_gpu
          REQUIRE_THROWS(from_fasta<uint8_t>(tmp_file_name, include_x,
                                             remove_duplicates, 2, use_gpu));
        } else {
          auto d = from_fasta<uint8_t>(tmp_file_name, include_x,
                                       remove_duplicates, 2, use_gpu);
          REQUIRE(d[{0, 0}] == 0);
          REQUIRE(d[{0, 1}] == 2);
          REQUIRE(d[{1, 0}] == 2);
          REQUIRE(d[{1, 1}] == 0);
        }
      }
    }
  }
  std::remove(tmp_file_name);
}

TEST_CASE("invalid input data: empty sequence", "[invalid]") {
  std::vector<std::vector<std::string>> expr;
  expr.push_back({});
  expr.push_back({{""}, {""}});
  std::string msg{"Error: Empty sequence"};
  for (auto &v : expr) {
    REQUIRE_THROWS_WITH(from_stringlist(v), msg);
  }
}

TEST_CASE("invalid input data: inconsistent sequence lengths", "[invalid]") {
  std::vector<std::vector<std::string>> expr;
  expr.push_back({{"ACGT"}, {"A"}});
  expr.push_back({{"AC"}, {"ACGTCG"}});
  expr.push_back({{"ACGT"}, {"ACGT"}, {"ACGT"}, {"ACGT"}, {"ACGT"}, {"A"}});
  expr.push_back({{"ACGT"}, {"A-----"}});
  expr.push_back({{"ACGT"}, {""}});
  std::string msg{"Error: Sequences do not all have the same length"};
  for (auto &v : expr) {
    REQUIRE_THROWS_WITH(from_stringlist(v), msg);
    REQUIRE_THROWS_WITH(distance(v.front(), v.back()), msg);
  }
}

TEST_CASE("from_csv reproduces correct data", "[hamming]") {
  std::mt19937 gen(12345);
  std::vector<std::string> data(107);
  for (auto &d : data)
    d = make_test_string(201, gen);

  DataSet<uint8_t> ref(data);
  char tmp_file_name[L_tmpnam];
  REQUIRE(std::tmpnam(tmp_file_name) != nullptr);
  ref.dump(std::string(tmp_file_name));

  auto restore = from_csv(std::string(tmp_file_name));
  REQUIRE(ref.nsamples == restore.nsamples);
  for (std::size_t i = 0; i < ref.nsamples; ++i) {
    for (std::size_t j = 0; j < ref.nsamples; ++j) {
      REQUIRE(ref[{i, j}] == restore[{i, j}]);
    }
  }
  std::remove(tmp_file_name);
}

TEMPLATE_TEST_CASE("distance integer saturates instead of overflowing",
                   "[hamming]", uint8_t, uint16_t) {
  auto n_max{static_cast<std::size_t>(std::numeric_limits<TestType>::max())};
  std::mt19937 gen(12345);
  std::vector<std::string> data(2);
  for (auto use_gpu : valid_use_gpu_values()) {
    for (auto n : {n_max, n_max + 1, n_max + 99}) {
      CAPTURE(use_gpu);
      CAPTURE(n);
      data[0] = std::string(n, 'A');
      data[1] = std::string(n, 'T');
      DataSet<TestType> dataSet(data, false, false, {}, use_gpu);
      REQUIRE(dataSet[{0, 1}] == n_max);
      REQUIRE(dataSet[{1, 0}] == n_max);
    }
  }
}

TEMPLATE_TEST_CASE("from_lower_triangular reproduces correct data", "[hamming]",
                   uint8_t, uint16_t) {
  std::mt19937 gen(12345);
  char tmp_file_name[L_tmpnam];
  REQUIRE(std::tmpnam(tmp_file_name) != nullptr);
  for (auto n : {1, 2, 3, 5, 9, 10, 19, 100, 207}) {
    CAPTURE(n);
    std::vector<std::string> data(n);
    for (auto &d : data)
      d = make_test_string(24, gen);

    DataSet<TestType> ref(data);
    REQUIRE(ref.nsamples == n);
    ref.dump_lower_triangular(std::string(tmp_file_name));

    auto restore = from_lower_triangular<TestType>(std::string(tmp_file_name));
    REQUIRE(ref.nsamples == restore.nsamples);
    for (std::size_t i = 0; i < ref.nsamples; ++i) {
      for (std::size_t j = 0; j < ref.nsamples; ++j) {
        REQUIRE(ref[{i, j}] == restore[{i, j}]);
      }
    }
  }
  std::remove(tmp_file_name);
}

TEST_CASE("from_stringlist GPU and CPU implementations give consistent results",
          "[hamming][gpu]") {
  if (!cuda_gpu_available()) {
    SKIP("No CUDA gpu available");
  }
  std::mt19937 gen(12345);
  for (bool include_x : {false}) {
    for (int n : {1, 5, 13, 32, 89, 185, 497, 1092}) {
      for (int n_samples : {2, 3, 4, 7, 11, 32, 33, 257, 689}) {
        std::vector<std::string> stringlist;
        stringlist.reserve(n_samples);
        for (std::size_t i = 0; i < n_samples; ++i) {
          stringlist.push_back(make_test_string(n, gen, include_x));
        }
        CAPTURE(include_x);
        auto d_cpu{from_stringlist(stringlist, include_x, false)};
        auto d_gpu{from_stringlist(stringlist, include_x, true)};
        for (std::size_t i = 0; i < n_samples; ++i) {
          for (std::size_t j = 0; j < n_samples; ++j) {
            REQUIRE(d_cpu[{i, j}] == d_gpu[{i, j}]);
          }
        }
      }
    }
  }
}

TEMPLATE_TEST_CASE(
    "from_fasta GPU and CPU implementations give consistent results",
    "[hamming][gpu]", uint8_t, uint16_t) {
  if (!cuda_gpu_available()) {
    SKIP("No CUDA gpu available");
  }
  char tmp_file_name[L_tmpnam];
  std::mt19937 gen(12345);
  REQUIRE(std::tmpnam(tmp_file_name) != nullptr);
  CAPTURE(tmp_file_name);
  for (bool remove_duplicates : {false, true}) {
    for (bool include_x : {false}) {
      for (int n : {17, 88, 381, 1023}) {
        for (int n_samples : {2, 3, 4, 7, 11, 127, 128, 255, 256, 257, 703}) {
          write_test_fasta(tmp_file_name, n, n_samples, gen, include_x);
          CAPTURE(include_x);
          auto d_cpu{from_fasta<TestType>(tmp_file_name, include_x,
                                          remove_duplicates, 0, false)};
          auto d_gpu{from_fasta<TestType>(tmp_file_name, include_x,
                                          remove_duplicates, 0, true)};
          for (std::size_t i = 0; i < n_samples; ++i) {
            for (std::size_t j = 0; j < n_samples; ++j) {
              REQUIRE(d_cpu[{i, j}] == d_gpu[{i, j}]);
            }
          }
        }
      }
    }
  }
  std::remove(tmp_file_name);
}

TEST_CASE("from_fasta_to_lower_triangular GPU consistent with CPU from_fasta",
          "[hamming][gpu]") {
  if (!cuda_gpu_available()) {
    SKIP("No CUDA gpu available");
  }
  std::mt19937 gen(12345);
  char tmp_fasta_file_name[L_tmpnam];
  REQUIRE(std::tmpnam(tmp_fasta_file_name) != nullptr);
  CAPTURE(tmp_fasta_file_name);
  char tmp_lt_file_name[L_tmpnam];
  REQUIRE(std::tmpnam(tmp_lt_file_name) != nullptr);
  CAPTURE(tmp_lt_file_name);
  for (bool remove_duplicates : {false, true}) {
    for (int n : {17, 88, 381, 1023}) {
      for (int n_samples :
           {2, 3, 4, 7, 11, 127, 128, 255, 256, 257, 703, 1012}) {
        CAPTURE(remove_duplicates);
        CAPTURE(n);
        CAPTURE(n_samples);
        write_test_fasta(tmp_fasta_file_name, n, n_samples, gen, false);
        auto d_cpu{from_fasta<uint16_t>(tmp_fasta_file_name, false,
                                        remove_duplicates, 0, false)};
        from_fasta_to_lower_triangular(tmp_fasta_file_name, tmp_lt_file_name,
                                       remove_duplicates, 0, true);
        auto d_gpu{from_lower_triangular<uint16_t>(tmp_lt_file_name)};
        for (std::size_t i = 0; i < n_samples; ++i) {
          for (std::size_t j = 0; j < n_samples; ++j) {
            CAPTURE(i);
            CAPTURE(j);
            REQUIRE(d_cpu[{i, j}] == d_gpu[{i, j}]);
          }
        }
      }
    }
  }
  std::remove(tmp_fasta_file_name);
  std::remove(tmp_lt_file_name);
}
