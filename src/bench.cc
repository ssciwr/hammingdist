#include "bench.hh"

#include <array>

namespace hamming {

std::string make_string(int64_t n, std::mt19937& gen){
    std::string s;
    s.reserve(n);
    std::uniform_int_distribution<> distrib(0, 4);
    std::array<char, 5> c{'A', 'C', 'G', 'T', '-'};
    for (int64_t i = 0; i < n; ++i) {
      s.push_back(c[distrib(gen)]);
    }
    return s;
}

std::vector<std::string> make_stringlist(int64_t n, std::mt19937& gen) {
  std::vector<std::string> v;
  v.reserve(n);
  for (int64_t row = 0; row < n; ++row) {
    v.push_back(make_string(n, gen));
  }
  return v;
}

}

BENCHMARK_MAIN();
