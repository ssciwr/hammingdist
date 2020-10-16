A small C++ tool to calculate pairwise distances between gene sequences given in fasta format.

# Prerequisites

The following software is currently required:

* git
* CMake >= 3.11
* A reasonably new C++ compiler

# Building

This sequence of commands lets you start from scratch:

```
git clone https://gitlab.dune-project.org/dominic/covid-tda.git
cd covid-tda
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make distance
```

This should (successfully) build the executable `distance` in the `build` subdirectory.
If `cmake` picks up the wrong compiler, it can be explicitly enforced by adding
`-DCMAKE_CXX_COMPILER=<path-to-compiler>` to the cmake call (in that case it is best
to remove the build directory and start from scratch).

# Running

The tool can be run from the command line:

```
./distance <path-to-input> <n>
```

Here, `<path-to-input>` must point to a fasta dataset, e.g. by putting `../data/example.fasta`
after putting a data file in the `data` directory. `n` is the maximum number of gene
sequences that the tool should read. This can be a smaller number than the number of
gene sequences in the dataset.

The output is currently written to a file `distances.csv`. The output is a full
matrix, not only the triangular part of the symmetric matrix.
