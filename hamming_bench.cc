#include "hamming.hh"
#include <benchmark/benchmark.h>
#include <random>
#include <string>
#include <vector>
#ifdef HAMMING_WITH_OPENMP
#include <omp.h>
#endif

static std::vector<std::string> make_stringlist(int64_t n) {
  std::vector<std::string> v;
  std::mt19937 gen(12345);
  std::uniform_int_distribution<> distrib(0, 4);
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
#ifdef HAMMING_WITH_OPENMP
  omp_set_num_threads(1);
#endif
  int64_t n{state.range(0)};
  auto v{make_stringlist(n)};
  for (auto _ : state) {
    from_stringlist(v);
  }
  state.SetComplexityN(n);
}

#ifdef HAMMING_WITH_OPENMP
static void bench_from_stringlist_omp(benchmark::State &state) {
  omp_set_num_threads(state.range(0));
  auto v{make_stringlist(4096)};
  for (auto _ : state) {
    from_stringlist(v);
  }
}
BENCHMARK(bench_from_stringlist_omp)->Arg(1)->Arg(2)->Arg(4)->Arg(8)->Arg(12)->Arg(24);
#endif

BENCHMARK(bench_from_stringlist)->RangeMultiplier(2)->Range(128, 4096)->Complexity();

BENCHMARK_MAIN();
