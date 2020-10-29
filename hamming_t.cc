#define CATCH_CONFIG_MAIN
#include "hamming.hh"
#include <catch2/catch.hpp>
#include <string>
#include <vector>

constexpr std::array<char, 4> valid_chars{'A', 'C', 'G', 'T'};

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

TEST_CASE("from_stringlist", "[hamming]") {
  std::vector<std::string> v{{"ACGT"}, {"AGGG"}};
  auto d = from_stringlist(v);
  REQUIRE(d[{0, 0}] == 0);
  REQUIRE(d[{0, 1}] == 2);
  REQUIRE(d[{1, 0}] == 2);
  REQUIRE(d[{1, 1}] == 0);
}
