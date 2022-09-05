#include "bench.hh"
#include "hamming/hamming.hh"
#include "hamming/hamming_impl.hh"
#ifdef HAMMING_WITH_OPENMP
#include <omp.h>
#endif

using namespace hamming;

static void bench_distance_cpp(benchmark::State &state) {
#ifdef HAMMING_WITH_OPENMP
  omp_set_num_threads(1);
#endif
  std::mt19937 gen(12345);
  int64_t n{state.range(0)};
  auto s1{from_string(make_string(n, gen))};
  auto s2{from_string(make_string(n, gen))};
  int d{0};
  for (auto _ : state) {
    d += distance_cpp(s1, s2);
  }
  state.SetComplexityN(n);
}

static void bench_distance_sparse(benchmark::State &state) {
#ifdef HAMMING_WITH_OPENMP
  omp_set_num_threads(1);
#endif
  std::mt19937 gen(12345);
  int64_t n{state.range(0)};
  auto s1{make_string(n, gen, false)};
  auto s2{s1};
  // make ~0.5% of s2 elements differ from s1
  randomize_n(s2, n / 200, gen);
  auto sparse = to_sparse_data({s1, s2}, false);
  int d{0};
  for (auto _ : state) {
    d += distance_sparse(sparse[0], sparse[1]);
  }
  state.SetComplexityN(n);
}

BENCHMARK(bench_distance_sparse)->Range(4096, 4194304)->Complexity();

BENCHMARK(bench_distance_cpp)->Range(4096, 4194304)->Complexity();
