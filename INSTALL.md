# Prerequisites

The following software is currently required if you build from scratch:

- git
- CMake >= 3.11
- A reasonably new C++ compiler
- Python 3
- OpenMP for shared memory parallelization (optional)

# Building

This sequence of commands lets you start from scratch:

```bash
git clone https://github.com/ssciwr/hammingdist.git
cd hammingdist
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make distance
make test
```

This should (successfully) build the executable `distance` in the `build` subdirectory,
and run the tests.
If `cmake` picks up the wrong compiler, it can be explicitly enforced by adding
`-DCMAKE_CXX_COMPILER=<path-to-compiler>` to the cmake call (in that case it is best
to remove the build directory and start from scratch).

To enable OpenMP multi-threading, add `-DHAMMING_WITH_OPENMP=ON` to the cmake call.

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

If running in parallel with OpenMP, you can control the number of threads available
by changing the environment variable `OMP_NUM_THREADS`:

```
OMP_NUM_THREADS=8 ./distance <path-to-input> <n>
```

# Building the Python interface

This sequence of command should build the Python interface:

```
git clone --recursive https://github.com/ssciwr/hammingdist.git
cd hammingdist

mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make hammingdist
```

If the wrong version of Python is found, it can be set by adding `-DPYTHON_EXECUTABLE=\path\to\python`
to the cmake command.

# Deploying the Python interface

The Python packages are currently built using Github Actions with the project
`ciwheelbuild`. A wheel build and deploy can be triggered by creating a release
in Github. Of course, Github Actions should be enabled and the PyPI API access
token needs to be stored a secret on the project.
