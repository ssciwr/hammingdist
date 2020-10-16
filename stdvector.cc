#include <vector>
#include <iostream>
#include <random>
#include <stdint.h>
#include <chrono>

using Gene = std::uint_fast8_t;

static std::vector<std::vector<Gene>> makeMatrix(std::size_t l, std::size_t n){
  std::vector<std::vector<Gene>> m;
  std::mt19937 gen(12345);
  std::uniform_int_distribution<> distrib(0, 5);
  m.reserve(n);
  for(std::size_t row=0; row<n; ++row){
    m.emplace_back();
    auto &v = m.back();
    v.reserve(l);
    for(std::size_t i=0; i<l; ++i){
      v.push_back(distrib(gen));
    }
  }
  return m;
}

static int distance(const std::vector<Gene>& a, const std::vector<Gene>& b){
  int r{0};
  for (std::size_t i=0; i<a.size(); ++i){
    if(a[i]!=b[i]){
      ++r;
    }
  }
  return r;
}

int main(int argc, char *argv[]	) {
  std::size_t length{30000};
  std::size_t nsamples = std::stoi(std::string(argv[1]));
  // std::string datapath = std::string(argv[2]);

  auto m = makeMatrix(length, nsamples);

  // Resize an output vector
  std::vector<std::vector<int>> result(nsamples);
  for(std::size_t i=0; i < nsamples; ++i)
    result[i].resize(i);

  // Call the distance function
  auto start = std::chrono::steady_clock::now();
  for(std::size_t i=0; i<nsamples; ++i)
    for(std::size_t j=0; j<i; ++j)
      result[i][j] = distance(m[i], m[j]);
  auto stop = std::chrono::steady_clock::now();

  // Write the results to the terminal
  for(std::size_t i=0; i < nsamples; ++i)
    for(std::size_t j=0; j<i; ++j)
      std::cout << result[i][j] << "\n";

  // Report measurements
  std::chrono::duration<double> elapsed_seconds = stop-start;
  std::cout << "Computation took: " << elapsed_seconds.count() << "s\n";

  return 0;
}
