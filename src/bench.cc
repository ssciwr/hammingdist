#include "bench.hh"

#include <array>
#include <fstream>

namespace hamming {

std::string make_string(int64_t n, std::mt19937 &gen, bool include_dash) {
  std::string s;
  s.reserve(n);
  std::size_t max{4};
  if (!include_dash) {
    max = 3;
  }
  std::uniform_int_distribution<> distrib(0, max);
  std::array<char, 5> c{'A', 'C', 'G', 'T', '-'};
  for (int64_t i = 0; i < n; ++i) {
    s.push_back(c[distrib(gen)]);
  }
  return s;
}

void randomize_n(std::string &str, std::size_t n, std::mt19937 &gen) {
  std::uniform_int_distribution<> distrib_loc(0, str.size() - 1);
  std::uniform_int_distribution<> distrib_char(0, 4);
  std::array<char, 5> c{'A', 'C', 'G', 'T', '-'};
  for (std::size_t i = 0; i < n; ++i) {
    str[distrib_loc(gen)] = c[distrib_char(gen)];
  }
}

std::vector<std::string> make_stringlist(int64_t n, int64_t n_samples,
                                         std::mt19937 &gen) {
  std::vector<std::string> v;
  v.reserve(n_samples);
  for (int64_t row = 0; row < n_samples; ++row) {
    v.push_back(make_string(n, gen));
  }
  return v;
}

void write_fasta(const std::string &filename, const std::string &seq,
                 std::size_t n_seq, std::mt19937 &gen,
                 std::size_t randomise_every_n) {
  std::ofstream fs;
  fs.open(filename);
  for (std::size_t i = 0; i < n_seq; ++i) {
    auto randomized_seq{seq};
    // randomly replace 1 in randomise_every_n elements of seq
    randomize_n(randomized_seq, seq.size() / randomise_every_n, gen);
    fs << ">seq" << i << "\n" << randomized_seq << "\n";
  }
  fs.close();
}

} // namespace hamming

BENCHMARK_MAIN();
