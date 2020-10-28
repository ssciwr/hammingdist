#include "hamming.hh"
#include <benchmark/benchmark.h>
#include <random>
#include <string>
#include <vector>

static std::vector<std::vector<Gene>> make_genes(int64_t n) {
  std::vector<std::vector<Gene>> v;
  std::mt19937 gen(12345);
  std::uniform_int_distribution<> distrib(0, 4);
  v.reserve(n);
  for (int64_t row = 0; row < n; ++row) {
    v.emplace_back();
    auto &r = v.back();
    r.reserve(n);
    for (int64_t i = 0; i < n; ++i) {
      r.push_back(distrib(gen));
    }
  }
  return v;
}

static std::vector<std::string>
to_stringlist(const std::vector<std::vector<Gene>> &genes) {
  std::vector<std::string> v;
  std::array<char, 4> c{'A', 'C', 'G', 'T'};
  for (const auto &gene : genes) {
    v.emplace_back();
    auto &r = v.back();
    for (const auto &g : gene) {
      r.push_back(c[g]);
    }
  }
  return v;
}

static void bench_from_stringlist(benchmark::State &state) {
  int64_t n{state.range(0)};
  auto v{to_stringlist(make_genes(n))};
  for (auto _ : state) {
    from_stringlist(v);
  }
  state.SetComplexityN(n);
}

static void bench_DataSet(benchmark::State &state) {
  int64_t n{state.range(0)};
  auto v{make_genes(n)};
  for (auto _ : state) {
    auto d = DataSet(v);
  }
  state.SetComplexityN(n);
}

BENCHMARK(bench_DataSet)->Range(32, 2048)->Complexity();
BENCHMARK(bench_from_stringlist)->Range(32, 2048)->Complexity();

BENCHMARK_MAIN();
