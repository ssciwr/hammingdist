A small C++ tool to calculate pairwise distances between gene sequences given in fasta format.

[![DOI](https://zenodo.org/badge/308676358.svg)](https://zenodo.org/badge/latestdoi/308676358)
[![pypi releases](https://img.shields.io/pypi/v/hammingdist.svg)](https://pypi.org/project/hammingdist)
[![python versions](https://img.shields.io/pypi/pyversions/hammingdist)](https://pypi.org/project/hammingdist)

# Python interface

To use the Python interface, you should install it from PyPI:

```bash
python -m pip install hammingdist
```

## Distances matrix

Then, you can e.g. use it in the following way from Python:

```python
import hammingdist

# To see the different optional arguments available:
help(hammingdist.from_fasta)

# To import all sequences from a fasta file
data = hammingdist.from_fasta("example.fasta")

# To import only the first 100 sequences from a fasta file
data = hammingdist.from_fasta("example.fasta", n=100)

# To import all sequences and remove any duplicates
data = hammingdist.from_fasta("example.fasta", remove_duplicates=True)

# To import all sequences from a fasta file, also treating 'X' as a valid character
data = hammingdist.from_fasta("example.fasta", include_x=True)

# The distance data can be accessed point-wise, though looping over all distances might be quite inefficient
print(data[14,42])

# The data can be written to disk in csv format (default `distance` Ripser format) and retrieved:
data.dump("backup.csv")
retrieval = hammingdist.from_csv("backup.csv")

# It can also be written in lower triangular format (comma-delimited row-major, `lower-distance` Ripser format):
data.dump_lower_triangular("lt.txt")
retrieval = hammingdist.from_lower_triangular("lt.txt")

# If the `remove_duplicates` option was used, the sequence indices can also be written.
# For each input sequence, this prints the corresponding index in the output:
data.dump_sequence_indices("indices.txt")

# Finally, we can pass the data as a list of strings in Python:
data = hammingdist.from_stringlist(["ACGTACGT", "ACGTAGGT", "ATTTACGT"])
```

## Duplicates

When `from_fasta` is called with the option `remove_duplicates=True`, duplicate sequences are removed before constructing the differences matrix.

For example given this set of three input sequences:

| Index | Sequence |
| ----- | -------- |
| 0     | ACG      |
| 1     | ACG      |
| 2     | TAG      |

The distances matrix would be a 2x2 matrix of distances between `ACG` and `TAG`:
| | ACG | TAG |
| --- | --- | --- |
| ACG | 0 | 2 |
| TAG | 2 | 0 |

The row of the distances matrix corresponding to each index in the original sequence would be:

| Index | Sequence | Row in distances matrix |
| ----- | -------- | ----------------------- |
| 0     | ACG      | 0                       |
| 1     | ACG      | 0                       |
| 2     | TAT      | 1                       |

This last column is what is written to disk by `DataSet.dump_sequence_indices`.

It can also be constructed (as a numpy array) without calculating the distances matrix by using `hammingdist.fasta_sequence_indices`

```python
import hammingdist

sequence_indices = hammingdist.fasta_sequence_indices(fasta_file)
```

## Large distance values

By default, the elements in the distances matrix returned by `hammingdist.from_fasta` have a maximum value of 255.

For distances larger than this `hammingdist.from_fasta_large` supports distances up to 65535 (but uses twice as much RAM)

## Distances from reference sequence

The distance of each sequence in a fasta file from a given reference sequence can be calculated using:

```python
import hammingdist

distances = hammingdist.fasta_reference_distances(sequence, fasta_file, include_x=True)
```

This function returns a numpy array that contains the distance of each sequence from the reference sequence.

You can also calculate the distance between two individual sequences:

```python
import hammingdist

distance = hammingdist.distance("ACGTX", "AAGTX", include_x=True)
```

## OpenMP on linux

The latest versions of hammingdist on linux are now built with OpenMP (multithreading) support.
If this causes any issues, you can install a previous version of hammingdist without OpenMP support:

```bash
pip install hammingdist==0.11.0
```
