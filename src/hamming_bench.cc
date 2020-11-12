#include "hamming/hamming.hh"
#include "hamming_impl.hh"
#include "bench.hh"
#ifdef HAMMING_WITH_OPENMP
#include <omp.h>
#endif

using namespace hamming;

static void bench_from_stringlist(benchmark::State &state) {
#ifdef HAMMING_WITH_OPENMP
  omp_set_num_threads(1);
#endif
  std::mt19937 gen(12345);
  int64_t n{state.range(0)};
  auto v{make_stringlist(n, gen)};
  for (auto _ : state) {
    from_stringlist(v);
  }
  state.SetComplexityN(n);
}

#ifdef HAMMING_WITH_OPENMP
static void bench_from_stringlist_omp(benchmark::State &state) {
  omp_set_num_threads(state.range(0));
  std::mt19937 gen(12345);
  auto v{make_stringlist(8192, gen)};
  for (auto _ : state) {
    from_stringlist(v);
  }
}
#endif

BENCHMARK(bench_from_stringlist)->RangeMultiplier(2)->Range(128, 8192)->Complexity();
#ifdef HAMMING_WITH_OPENMP
BENCHMARK(bench_from_stringlist_omp)->Arg(1)->Arg(2)->Arg(4)->Arg(8)->Arg(12)->Arg(24);
#endif
