#include <array>
#include <vector>
#include <iostream>
#include <stdint.h>
#include <chrono>
#include <fstream>
#include <sstream>

using Gene = std::uint_fast8_t;

static std::vector<std::vector<Gene>> makeMatrix(const std::string& filename, std::size_t l, std::size_t n)
{
  // The return variable
  std::vector<std::vector<Gene>> m;
  m.reserve(n);

  // Initializing the stream
  std::ifstream stream(filename);
  std::size_t count = 0;
  std::string line;
  unsigned char buffer;

  // Lookup table for the encodings
  std::array<Gene, 256> lookup;
  std::fill(lookup.begin(), lookup.end(), 0);
  lookup[std::size_t('A')] = 1 << 1;
  lookup[std::size_t('C')] = 1 << 2;
  lookup[std::size_t('G')] = 1 << 3;
  lookup[std::size_t('T')] = 1 << 4;

  // The while-condition consumes the header line
  while((count < n) && (std::getline(stream, line)))
  {
    m.emplace_back();
    auto &v = m.back();
    v.resize(l);
    std::size_t pos = 0;

    while((pos < l) && (std::getline(stream, line)))
    {
      std::istringstream iss(line);
      while(iss >> buffer)
        v[pos++] = lookup[std::size_t(buffer)];
    }
   
    ++count;
  }

  return m;
}

static int distance(const std::vector<Gene>& a, const std::vector<Gene>& b){
  int r{0};
  for (std::size_t i=0; i<a.size(); ++i)
    if(a[i] ^ b[i])
      ++r;
  return r;
}

int main(int argc, char *argv[]	) {
  // Parse arguments. Sequence length currently fixed at compile time
  std::size_t length{29768};
  std::string filename = std::string(argv[1]);
  std::size_t nsamples = std::stoi(std::string(argv[2]));

  // Read the input
  auto start = std::chrono::steady_clock::now();
  auto m = makeMatrix(filename, length, nsamples);
  auto stop = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = stop - start;
  std::cout << "Parsing input fasta file took " << elapsed_seconds.count() << "s" << std::endl;

  // Properly resize the output array
  start = std::chrono::steady_clock::now();
  std::vector<std::vector<int>> result(nsamples);
  for(std::size_t i=0; i < nsamples; ++i)
    result[i].resize(i);
  stop = std::chrono::steady_clock::now();
  elapsed_seconds = stop - start;
  std::cout << "Allocating the output data structure took " << elapsed_seconds.count() << "s" << std::endl;

  // Call the distance function
  start = std::chrono::steady_clock::now();
  for(std::size_t i=0; i<nsamples; ++i)
    for(std::size_t j=0; j<i; ++j)
      result[i][j] = distance(m[i], m[j]);
  stop = std::chrono::steady_clock::now();
  elapsed_seconds = stop - start;
  std::cout << "Calculating distances took " << elapsed_seconds.count() << "s" << std::endl;

  // Write output to a file
  start = std::chrono::steady_clock::now();
  std::string outputfile{"distances.csv"};
  std::ofstream stream(outputfile);
  for(std::size_t i=0; i<nsamples; ++i)
  {
    for(std::size_t j=0; j<nsamples; ++j)
    {
      if(i == j)
        stream << 0;
      if(i < j)
        stream << result[j][i];
      if(i > j)
        stream << result[i][j];
      if (j != nsamples - 1)
        stream << ", ";
    }
    stream << std::endl;
  }
  stop = std::chrono::steady_clock::now();
  elapsed_seconds = stop - start;
  std::cout << "Writing results to " << outputfile << " took: " << elapsed_seconds.count() << "s" << std::endl;

  return 0;
}
 