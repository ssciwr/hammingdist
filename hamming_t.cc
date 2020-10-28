#define CATCH_CONFIG_MAIN
#include "hamming.hh"
#include <catch2/catch.hpp>
#include <string>
#include <vector>

TEST_CASE("from_stringlist", "[hamming]") {
  std::vector<std::string> v{{"ACGT"}, {"AGGG"}};
  auto d = from_stringlist(v);
  REQUIRE(d[{0, 0}] == 0);
  REQUIRE(d[{0, 1}] == 2);
  REQUIRE(d[{1, 0}] == 2);
  REQUIRE(d[{1, 1}] == 0);
}
