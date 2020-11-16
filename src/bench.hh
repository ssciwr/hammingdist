#ifndef HAMMING_BENCH_HH
#define HAMMING_BENCH_HH

#include <benchmark/benchmark.h>
#include <string>
#include <random>
#include <vector>

namespace hamming {

std::string make_string(int64_t n, std::mt19937& gen, bool include_dash = true);
void randomize_n(std::string& str, std::size_t n, std::mt19937& gen);
std::vector<std::string> make_stringlist(int64_t n, std::mt19937& gen);

}

#endif
