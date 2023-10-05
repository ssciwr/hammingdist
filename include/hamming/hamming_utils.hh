#pragma once

#include <fmt/core.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#ifdef HAMMING_WITH_OPENMP
#include <omp.h>
#endif

namespace hamming {

template <typename IntDistType>
void append_int(std::string &str, IntDistType value) {
  str.append(fmt::to_string(value));
}

template <typename DistIntType>
std::string
lower_triangular_lines(const std::vector<DistIntType> &partial_distances,
                       std::size_t i_start, std::size_t i_end, std::size_t i0,
                       std::size_t j0, std::size_t iN, std::size_t jN) {
  std::string lines{};
  std::size_t k{[i_start, i0, j0]() {
    if (i_start == i0) {
      // first row, no offset required
      return std::size_t{0};
    }
    // offset for elements from first row and any previous full rows
    return i0 - j0 + i_start * (i_start - 1) / 2 - (i0 + 1) * i0 / 2;
  }()};
  for (std::size_t i = i_start; i < i_end; ++i) {
    // first row starts at j0, other rows are full and start at 0
    std::size_t j_first{i == i0 ? j0 : 0};
    std::size_t j_last{i - 1};
    char final_char{'\n'};
    if (i == iN && j_last != jN) {
      // this is an incomplete final row of partial_distances: no newline at end
      j_last = jN;
      final_char = ',';
    }
    // do all but last element for this line
    for (std::size_t j = j_first; j < j_last; ++j) {
      append_int(lines, partial_distances[k++]);
      lines.push_back(',');
    }
    // do last element for this line with ',' or '\n' as appropriate
    append_int(lines, partial_distances[k++]);
    lines.push_back(final_char);
  }
  return lines;
}

static std::size_t row_from_index(std::size_t index) {
  return static_cast<std::size_t>(
      floor(sqrt(2.0 * static_cast<double>(index) + 0.5) + 0.5));
}

static std::size_t col_from_index(std::size_t index, std::size_t row) {
  return index - row * (row - 1) / 2;
}

template <typename DistIntType>
void partial_write_lower_triangular(
    const std::string &filename,
    const std::vector<DistIntType> &partial_distances,
    std::size_t distances_offset, std::size_t n_partial_distances) {
  if (n_partial_distances == 0) {
    return;
  }
  // infer row/col (i0, j0) of first element of partial_distances
  std::size_t i0{row_from_index(distances_offset)};
  std::size_t j0{col_from_index(distances_offset, i0)};
  // infer row/col (iN, jN) of last element of partial_distances
  std::size_t last_element_index{distances_offset + n_partial_distances - 1};
  std::size_t iN{row_from_index(last_element_index)};
  std::size_t jN{col_from_index(last_element_index, iN)};
  auto flag = [distances_offset]() {
    if (distances_offset == 0) {
      // overwrite any existing data
      return std::ios_base::out;
    }
    // append to existing data
    return std::ios_base::out | std::ios_base::app | std::ios_base::ate;
  }();
  std::ofstream output_file_stream(filename, flag);
#ifdef HAMMING_WITH_OPENMP
  constexpr std::size_t max_lines_per_thread{200};
#pragma omp parallel for schedule(static, 1) ordered default(none)             \
    shared(partial_distances, output_file_stream, i0, j0, iN, jN)
  for (std::size_t i_start = i0; i_start <= iN;
       i_start += max_lines_per_thread) {
    std::size_t i_end{std::min(i_start + max_lines_per_thread, iN + 1)};
    auto line = lower_triangular_lines(partial_distances, i_start, i_end, i0,
                                       j0, iN, jN);
#pragma omp ordered
    output_file_stream << line;
  }
#else
  output_file_stream << lower_triangular_lines(partial_distances, i0, iN + 1,
                                               i0, j0, iN, jN);
#endif
}

template <typename DistIntType>
void write_lower_triangular(const std::string &filename,
                            const std::vector<DistIntType> &distances) {
  partial_write_lower_triangular(filename, distances, 0, distances.size());
}

} // namespace hamming
