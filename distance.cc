#include"hamming/hamming.hh"

#include<chrono>
#include<iostream>
#include<string>

int main(int argc, char *argv[]	) {
  std::string filename = std::string(argv[1]);
  std::size_t nsamples = std::stoi(std::string(argv[2]));

  auto start = std::chrono::steady_clock::now();
  auto data = hamming::from_fasta(filename, nsamples);
  data.dump("distances.csv");
  auto stop = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = stop - start;
  std::cout << "Processing fasta data set took " << elapsed_seconds.count() << "s" << std::endl;
  return 0;
}
 
