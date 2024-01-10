#include "cuda_mem.hh"
#include "hamming/distance_cuda.hh"
#include "hamming/hamming_utils.hh"
#include <chrono>
#include <cuda/std/limits>

namespace hamming {

template <typename DistIntType>
__global__ void Dist(DistIntType *partial_distances, const std::uint8_t *genes,
                     std::uint64_t distances_offset,
                     unsigned int geneBlocksPerSample) {
  // Calculates all gridDim.x entries of the partial_distances array.
  //
  // The full distances array is a flat nsamples * (nsamples - 1) / 2 element
  // array that contains the lower-triangular elements of the nSamples x
  // nSamples partial_distances matrix.
  //
  // This kernel is provided with partial_distances, which should have
  // gridDim.x entries, and which is filled with values corresponding to
  // entries in the full distances array with an offset of distances_offset.
  //

  // this array is shared between the threads in this block
  // and must be large enough to store one int per thread
  extern __shared__ int s[];

  // index in the partial_distances array where we'll put the result from this
  // block of threads
  uint64_t distancesIndex{static_cast<uint64_t>(blockIdx.x)};
  // index of this value in the full distances array
  uint64_t trueDistancesIndex{distancesIndex + distances_offset};
  // infer indices of the two genes corresponding to this distances index
  uint64_t distancesRowIndex{static_cast<std::size_t>(
      floor(sqrt(2.0 * static_cast<double>(trueDistancesIndex) + 0.5) + 0.5))};
  uint64_t distancesColIndex{trueDistancesIndex -
                             distancesRowIndex * (distancesRowIndex - 1) / 2};
  uint64_t uint32sPerSample{geneBlocksPerSample / 4};
  uint64_t geneAIndex{distancesRowIndex * uint32sPerSample};
  uint64_t geneBIndex{distancesColIndex * uint32sPerSample};

  unsigned int threadIndex{threadIdx.x};
  // calculate partial sum for each thread and store in shared memory s
  int r0{0};
  int r1{0};
  int r2{0};
  int r3{0};
  // NOTE: this cast is only safe if genes is 32-bit aligned AND we each sample
  // in genes is padded such that the first element of each sample is also
  // 32-bit aligned!
  const uint32_t *genes_as_uint32{reinterpret_cast<const uint32_t *>(genes)};
  uint mask_lower{0x0f0f0f0f};
  uint mask_upper{0xf0f0f0f0};
  // NOTE: this loop is also only correct if the length of each sample in genes
  // is a multiple of 8, which we do by padding the samples with '-'.
  for (int j = 2 * threadIndex; j < uint32sPerSample; j += 2 * blockDim.x) {
    auto c0{genes_as_uint32[geneAIndex + j] & genes_as_uint32[geneBIndex + j]};
    auto c1{genes_as_uint32[geneAIndex + j + 1] &
            genes_as_uint32[geneBIndex + j + 1]};
    r0 += __popc(__vseteq4(c0 & mask_lower, 0u));
    r1 += __popc(__vseteq4(c0 & mask_upper, 0u));
    r2 += __popc(__vseteq4(c1 & mask_lower, 0u));
    r3 += __popc(__vseteq4(c1 & mask_upper, 0u));
  }
  s[threadIndex] = r0 + r1 + r2 + r3;
  // synchronise shared memory s between all threads in this block
  __syncthreads();
  // sum elements of s using reduction until partial sums are stored in
  // the first 64 elements of s
  for (int offset = blockDim.x / 2; offset > 32; offset >>= 1) {
    if (threadIndex < offset) {
      s[threadIndex] += s[threadIndex + offset];
    }
    __syncthreads();
  }
  if (threadIndex < 32) {
    // one more reduction in each of the 32 threads in this warp
    int sum{s[threadIndex] + s[threadIndex + 32]};
    // now sum the values of sum within this warp
    constexpr unsigned int FULL_MASK{0xffffffff};
    for (int offset = 16; offset > 0; offset /= 2) {
      sum += __shfl_down_sync(FULL_MASK, sum, offset);
    }
    if (threadIndex == 0) {
      auto maxDist{cuda::std::numeric_limits<DistIntType>::max()};
      partial_distances[distancesIndex] = sum > maxDist ? maxDist : sum;
    }
  }
}

template <typename DistIntType>
std::vector<DistIntType>
distances_cuda(const std::vector<std::vector<GeneBlock>> &data,
               const std::string &filename = {}) {
  std::vector<DistIntType> distances{};
  std::size_t timing_gpu_ms = 0;
  std::size_t timing_io_ms = 0;
  auto timing0{std::chrono::high_resolution_clock::now()};
  bool output_to_vector{false};
  if (filename.empty()) {
    output_to_vector = true;
  }
  std::size_t nSamples{data.size()};
  std::size_t nDistances{nSamples * (nSamples - 1) / 2};
  std::size_t geneBlocksPerSample{data[0].size()};
  // 2^31-1 is limit on number of CUDA blocks in x-dim, which corresponds to
  // 0.5/1GB of distances data for each chunk. For large datasets I/O becomes
  // the bottleneck and the larger the chunk the faster the I/O tends to be.
  std::size_t nPartialDistances{std::min(nDistances, 2147483647ul)};

  if (output_to_vector) {
    // need to store all distances on host
    distances.resize(nDistances);
  } else {
    // only need to store a single block of partial distances on host
    distances.resize(nPartialDistances);
  }
  // allocate memory for genes on device
  // one gene is 30k chars -> 15k bytes in dense format
  // so 1 million samples -> 15GB
  auto *genes{CheckedCudaMalloc<GeneBlock>(nSamples * geneBlocksPerSample)};
  // copy genes to device
  for (std::size_t i = 0; i < data.size(); ++i) {
    CheckedCopyToDevice(genes + i * geneBlocksPerSample, data[i]);
  }

  // allocate memory for partial distances matrix on device
  auto *partial_distances{CheckedCudaMalloc<DistIntType>(nPartialDistances)};
  // keep track of how many distance elements are available to write to disk
  std::size_t available_distance_elements{0};
  // keep track of where in the full distances array these elements should go
  std::size_t distances_offset{0};
  // GPU timing
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  // I/O timing
  auto timing_io = std::chrono::high_resolution_clock::now();
  while (distances_offset + available_distance_elements < nDistances) {
    // use nThreadsPerBlock in x dim of block
    uint nThreadsPerBlock{128};
    dim3 threadsPerBlock{nThreadsPerBlock, 1, 1};
    // use up to nPartialDistances blocks, one block per distance element
    dim3 numBlocks{static_cast<uint>(std::min(nPartialDistances,
                                              nDistances - distances_offset)),
                   1, 1};
    cudaEventRecord(start);
    timing_io = std::chrono::high_resolution_clock::now();
    // launch a kernel with shared memory of size int[nThreadsPerBlock] -
    // this call returns immediately and the kernel runs asynchronously on the
    // GPU
    Dist<<<numBlocks, threadsPerBlock, nThreadsPerBlock * sizeof(int)>>>(
        partial_distances, genes, distances_offset, geneBlocksPerSample);
    cudaEventRecord(stop);
    if (auto err = cudaGetLastError(); err != cudaSuccess) {
      throw std::runtime_error(cudaGetErrorString(err));
    }

    if (!output_to_vector) {
      // write previous kernel's output (if any) to disk using CPU while the new
      // kernel is running on the GPU to interleave I/O with computation
      partial_write_lower_triangular(filename, distances, distances_offset,
                                     available_distance_elements);
    }
    timing_io_ms += std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::high_resolution_clock::now() - timing_io)
                        .count();
    // copy partial_distances from GPU to distances vector on HOST - this call
    // waits until the kernel has completed before copying the memory.
    CheckedCopyToHost(distances.data() +
                          (output_to_vector ? distances_offset : 0),
                      partial_distances, numBlocks.x);
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    timing_gpu_ms += static_cast<std::size_t>(milliseconds);
    distances_offset += available_distance_elements;
    available_distance_elements = numBlocks.x;
  }
  if (!output_to_vector) {
    // write final kernel's output to disk
    timing_io = std::chrono::high_resolution_clock::now();
    partial_write_lower_triangular(filename, distances, distances_offset,
                                   available_distance_elements);
    timing_io_ms += std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::high_resolution_clock::now() - timing_io)
                        .count();
    distances.clear();
  }
  // free data on gpu
  cudaFree(genes);
  cudaFree(partial_distances);
  std::cout << "# hammingdist :: ...distance calculation completed in "
            << std::chrono::duration_cast<std::chrono::milliseconds>(
                   std::chrono::high_resolution_clock::now() - timing0)
                   .count()
            << " ms (GPU: " << timing_gpu_ms << " / IO: " << timing_io_ms
            << ")." << std::endl;
  return distances;
}

std::vector<uint8_t>
distances_cuda_8bit(const std::vector<std::vector<GeneBlock>> &data) {
  return distances_cuda<uint8_t>(data, {});
}

std::vector<uint16_t>
distances_cuda_16bit(const std::vector<std::vector<GeneBlock>> &data) {
  return distances_cuda<uint16_t>(data, {});
}

void distances_cuda_to_lower_triangular(
    const std::vector<std::vector<GeneBlock>> &data,
    const std::string &filename) {
  distances_cuda<uint16_t>(data, filename);
}

int distance_cuda(const std::vector<GeneBlock> &a,
                  const std::vector<GeneBlock> &b) {
  // wrapper for testing cuda kernel with existing distance API
  std::vector<std::vector<GeneBlock>> data{a, b};
  return distances_cuda<int>(data, {})[0];
}

bool distance_cuda_have_device() {
  int nDevices = 0;
  if (cudaError_t err{cudaGetDeviceCount(&nDevices)}; err != cudaSuccess) {
    return false;
  }
  return nDevices > 0;
}

} // namespace hamming
