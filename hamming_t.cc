#define CATCH_CONFIG_MAIN
#include "hamming.hh"
#include <catch2/catch.hpp>
#include <string>
#include <vector>
#include <fstream>
#include <cstdio>

constexpr std::array<char, 4> valid_chars{'A', 'C', 'G', 'T'};
constexpr std::array<char, 6> invalid_chars{' ', 'N', '*', '?', 'a', '.'};

static int dist(char c1, char c2) {
  return from_stringlist({std::string{c1}, std::string{c2}})[{0, 1}];
}

TEST_CASE("distance between two equal valid characters is 0", "[distance]") {
  for (auto c : valid_chars) {
    CAPTURE(c);
    REQUIRE(dist(c, c) == 0);
  }
}

TEST_CASE("distance between any valid char & '-' is 0", "[distance]") {
  for (auto c : valid_chars) {
    CAPTURE(c);
    REQUIRE(dist(c, '-') == 0);
    REQUIRE(dist('-', c) == 0);
  }
}

TEST_CASE("distance between two different valid characters is 1",
          "[distance]") {
  for (auto c1 : valid_chars) {
    for (auto c2 : valid_chars) {
      if (c1 != c2) {
        CAPTURE(c1);
        CAPTURE(c2);
        REQUIRE(dist(c1, c2) == 1);
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
      REQUIRE(dist(c1, c2) == 1);
    }
  }
}

TEST_CASE("distance between two invalid characters is 1", "[distance]") {
  for (auto c1 : invalid_chars) {
    for (auto c2 : invalid_chars) {
      CAPTURE(c1);
      CAPTURE(c2);
      REQUIRE(dist(c1, c2) == 1);
    }
  }
}

TEST_CASE("distance between any invalid char & '-' is 1", "[distance]") {
  for (auto c : invalid_chars) {
    CAPTURE(c);
    REQUIRE(dist(c, '-') == 1);
    REQUIRE(dist('-', c) == 1);
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
  expr.push_back({{"ACGTGTCGTGTCGACGTGTCG---"}, {"ACGTGTCGTTTCGACGAGTCGGGG"}});
  for (const auto &v : expr) {
    auto d = from_stringlist(v);
    REQUIRE(d[{0, 0}] == 0);
    REQUIRE(d[{0, 1}] == 2);
    REQUIRE(d[{1, 0}] == 2);
    REQUIRE(d[{1, 1}] == 0);
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
    auto d = from_fasta(tmp_file_name, 2);
    REQUIRE(d[{0, 0}] == 0);
    REQUIRE(d[{0, 1}] == 2);
    REQUIRE(d[{1, 0}] == 2);
    REQUIRE(d[{1, 1}] == 0);
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
    auto d = from_fasta(tmp_file_name, 2);
    REQUIRE(d[{0, 0}] == 0);
    REQUIRE(d[{0, 1}] == 2);
    REQUIRE(d[{1, 0}] == 2);
    REQUIRE(d[{1, 1}] == 0);
    std::remove(tmp_file_name);
}

TEST_CASE("invalid input data: empty sequence", "[invalid]") {
    std::vector<std::vector<std::string>> expr;
    expr.push_back({});
    expr.push_back({{""}, {""}});
    std::string msg{"Error: Empty sequence"};
    for (const auto &v : expr) {
        REQUIRE_THROWS_WITH(from_stringlist(v), msg);
    }
}

TEST_CASE("invalid input data: inconsistent sequence lengths", "[invalid]") {
    std::vector<std::vector<std::string>> expr;
    expr.push_back({{"ACGT"}, {"A"}});
    expr.push_back({{"AC"}, {"ACGTCG"}});
    expr.push_back({{"ACGT"}, {"ACGT"}, {"ACGT"}, {"A"}, {"ACGT"}, {"ACGT"}});
    expr.push_back({{"ACGT"}, {"A-----"}});
    expr.push_back({{"ACGT"}, {""}});
    std::string msg{"Error: Sequences do not all have the same length"};
    for (const auto &v : expr) {
        REQUIRE_THROWS_WITH(from_stringlist(v), msg);
    }
}
