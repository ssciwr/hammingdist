A small C++ tool to calculate pairwise distances between gene sequences given in fasta format.

[![DOI](https://zenodo.org/badge/308676358.svg)](https://zenodo.org/badge/latestdoi/308676358)
[![pypi releases](https://img.shields.io/pypi/v/hammingdist.svg)](https://pypi.org/project/hammingdist)
[![python versions](https://img.shields.io/pypi/pyversions/hammingdist)](https://pypi.org/project/hammingdist)

# Python interface

To use the Python interface, you should install it from PyPI:

```bash
python -m pip install hammingdist
```

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

## Distances from reference sequence

The distance of each sequence in a fasta file from a given reference sequence can be calculated using:

```python
import hammingdist

distances = hammingdist.fasta_reference_distances(sequence, fasta_file, include_x=True)
```

This function returns a numpy array that contains the distance of each sequence from the reference sequence.

## OpenMP on linux

The latest version of hammingdist on linux is now built with OpenMP (multithreading) support.
If this causes any issues, you can install a previous version of hammingdist without OpenMP support:

```bash
pip install hammingdist==0.11.0
```
