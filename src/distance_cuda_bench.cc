#include "bench.hh"
#include "hamming/distance_cuda.hh"
#include "hamming/hamming.hh"

using namespace hamming;

static void bench_distance_cuda(benchmark::State &state) {
  std::mt19937 gen(12345);
  int64_t n{state.range(0)};
  auto s1{from_string(make_string(n, gen))};
  auto s2{from_string(make_string(n, gen))};
  int d{0};
  for (auto _ : state) {
    d += distance_cuda(s1, s2);
  }
  state.SetComplexityN(n);
}

BENCHMARK(bench_distance_cuda)->Range(4096, 4194304)->Complexity();
