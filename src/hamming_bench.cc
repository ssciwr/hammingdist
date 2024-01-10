#include "bench.hh"
#include "hamming/hamming.hh"
#include "hamming/hamming_impl.hh"
#ifdef HAMMING_WITH_OPENMP
#include <omp.h>
#endif

using namespace hamming;

constexpr int64_t sampleLength{30000};

static void bench_from_stringlist(benchmark::State &state) {
#ifdef HAMMING_WITH_OPENMP
  omp_set_num_threads(1);
#endif
  std::mt19937 gen(12345);
  int64_t n{state.range(0)};
  auto v{make_stringlist(sampleLength, n, gen)};
  for (auto _ : state) {
    from_stringlist(v);
  }
  state.SetComplexityN(n);
}

#ifdef HAMMING_WITH_OPENMP
static void bench_from_stringlist_omp(benchmark::State &state) {
  omp_set_num_threads(state.range(0));
  std::mt19937 gen(12345);
  auto v{make_stringlist(sampleLength, 1024, gen)};
  for (auto _ : state) {
    from_stringlist(v);
  }
}
#endif

static void bench_from_stringlist_gpu(benchmark::State &state) {
  std::mt19937 gen(12345);
  int64_t n{state.range(0)};
  auto v{make_stringlist(sampleLength, n, gen)};
  for (auto _ : state) {
    from_stringlist(v, false, true);
  }
  state.SetComplexityN(n);
}

static void bench_from_fasta_to_lower_triangular_gpu(benchmark::State &state) {
  std::mt19937 gen(12345);
  int64_t n{state.range(0)};
  std::string fasta_file{benchmark_tmp_input_file};
  std::string lt_file{benchmark_tmp_output_file};
  auto reference_seq{make_string(sampleLength, gen, true)};
  write_fasta(fasta_file, reference_seq, state.range(0), gen);
  for (auto _ : state) {
    from_fasta_to_lower_triangular(fasta_file, lt_file, false);
  }
  state.SetComplexityN(n);
}

static void bench_fasta_reference_distances(benchmark::State &state) {
  std::mt19937 gen(12345);
  std::string fasta_file{benchmark_tmp_input_file};
  auto reference_seq{make_string(sampleLength, gen, true)};
  write_fasta(fasta_file, reference_seq, state.range(0), gen);
  std::vector<ReferenceDistIntType> distances;
  for (auto _ : state) {
    distances = fasta_reference_distances(reference_seq, fasta_file, true);
  }
}

BENCHMARK(bench_from_stringlist)
    ->RangeMultiplier(2)
    ->Range(16, 1024)
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
#ifdef HAMMING_WITH_CUDA
BENCHMARK(bench_from_stringlist_gpu)
    ->RangeMultiplier(2)
    ->Range(16, 8192)
    ->Complexity();
BENCHMARK(bench_from_fasta_to_lower_triangular_gpu)
    ->RangeMultiplier(2)
    ->Range(16, 16384)
    ->Complexity();
#endif
BENCHMARK(bench_fasta_reference_distances)
    ->RangeMultiplier(2)
    ->Range(16, 1024)
    ->Complexity();
