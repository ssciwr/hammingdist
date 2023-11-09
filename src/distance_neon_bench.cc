#include "bench.hh"
#include "hamming/distance_neon.hh"
#include "hamming/hamming.hh"
#include "hamming/hamming_impl.hh"
#ifdef HAMMING_WITH_OPENMP
#include <omp.h>
#endif

using namespace hamming;

static void bench_distance_neon(benchmark::State &state) {
#ifdef HAMMING_WITH_OPENMP
  omp_set_num_threads(1);
#endif
  std::mt19937 gen(12345);
  int64_t n{state.range(0)};
  auto s1{from_string(make_string(n, gen))};
  auto s2{from_string(make_string(n, gen))};
  int d{0};
  for (auto _ : state) {
    d += distance_neon(s1, s2);
  }
  state.SetComplexityN(n);
}

BENCHMARK(bench_distance_neon)->Range(4096, 4194304)->Complexity();
