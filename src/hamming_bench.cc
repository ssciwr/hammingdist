#include "bench.hh"
#include "hamming/hamming.hh"
#include "hamming/hamming_impl.hh"
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

static void bench_fasta_reference_distances(benchmark::State &state) {
  std::mt19937 gen(12345);
  std::string fasta_file{"fasta.txt"};
  auto reference_seq{make_string(30000, gen, true)};
  write_fasta(fasta_file, reference_seq, state.range(0), gen);
  std::vector<ReferenceDistIntType> distances;
  for (auto _ : state) {
    distances = fasta_reference_distances(reference_seq, fasta_file, true);
  }
}

BENCHMARK(bench_from_stringlist)
    ->RangeMultiplier(2)
    ->Range(128, 8192)
    ->Complexity();
#ifdef HAMMING_WITH_OPENMP
BENCHMARK(bench_from_stringlist_omp)
    ->Arg(1)
    ->Arg(2)
    ->Arg(4)
    ->Arg(8)
    ->Arg(12)
    ->Arg(24);
#endif
BENCHMARK(bench_fasta_reference_distances)
    ->RangeMultiplier(2)
    ->Range(16, 32384)
    ->Complexity();
