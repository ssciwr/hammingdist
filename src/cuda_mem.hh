#pragma once

#include <stdexcept>
#include <vector>

template <typename T> T *CheckedCudaMalloc(std::size_t n) {
  T *a{nullptr};
  std::size_t sz{sizeof(T) * n};
  if (cudaError err{cudaMalloc(&a, sz)}; err != cudaSuccess) {
    throw std::runtime_error(cudaGetErrorString(err));
  }
  return a;
}

template <typename T>
void CheckedCopyToDevice(T *dest, const std::vector<T> &src) {
  std::size_t count{sizeof(T) * src.size()}; // count in bytes
  if (cudaError err{
          cudaMemcpy(dest, src.data(), count, cudaMemcpyHostToDevice)};
      err != cudaSuccess) {
    throw std::runtime_error(cudaGetErrorString(err));
  }
}

template <typename T>
void CheckedCopyToHost(T *dest, const T *src, std::size_t n_elements) {
  std::size_t count{sizeof(T) * n_elements}; // count in bytes
  if (cudaError err{cudaMemcpy(dest, src, count, cudaMemcpyDeviceToHost)};
      err != cudaSuccess) {
    throw std::runtime_error(cudaGetErrorString(err));
  }
}
