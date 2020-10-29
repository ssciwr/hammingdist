#include "hamming.hh"
#include <benchmark/benchmark.h>
#include <random>
#include <string>
#include <vector>

static std::vector<std::string> make_stringlist(int64_t n) {
  std::vector<std::string> v;
  std::mt19937 gen(12345);
  std::uniform_int_distribution<> distrib(0, 5);
  std::array<char, 5> c{'A', 'C', 'G', 'T', '-'};
  v.reserve(n);
  for (int64_t row = 0; row < n; ++row) {
    v.emplace_back();
    auto &r = v.back();
    r.reserve(n);
    for (int64_t i = 0; i < n; ++i) {
      r.push_back(c[distrib(gen)]);
    }
  }
  return v;
}

static void bench_from_stringlist(benchmark::State &state) {
  int64_t n{state.range(0)};
  auto v{make_stringlist(n)};
  for (auto _ : state) {
    from_stringlist(v);
  }
  state.SetComplexityN(n);
}

BENCHMARK(bench_from_stringlist)->Range(32, 2048)->Complexity();

BENCHMARK_MAIN();
