A small C++ tool to calculate pairwise distances between gene sequences given in fasta format.

# Python interface

To use the Python interface, you should install it from PyPI:

```
python -m pip install hammingdist
```

Then, you can e.g. use it in the following way from Python:

```
import hammingdist

# This accepts excactly the same two arguments as the stand-alone
# executable: A fasta file and the maximum number of sequences to consider
data = hammingdist.from_fasta("example.fasta", 100)

# The distance data can be accessed point-wise, though looping over all distances might be quite inefficient
print(data[14,42])

# The data can be written to disk and retrieved:
data.dump("backup.csv")
retrieval = hammingdist.from_csv("backup.csv")

# Finally, we can pass the data as a list of strings in Python:
data = hammingdist.from_stringlist(["ACGTACGT", "ACGTAGGT", "ATTTACGT"])

# When in doubt, the internal data structures of the DataSet object can be inspected:
print(data._data)
print(data._distances)
```

# Prerequisites

The following software is currently required if you build from scratch:

* git
* CMake >= 3.11
* A reasonably new C++ compiler
* Python 3

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

# Building the Python interface

This sequence of command should build the Python interface:

```
git clone --recursive https://gitlab.dune-project.org/dominic/covid-tda.git
cd covid-tda

mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make hammingdist
```

# Deploying the Python interface

In order to do this, `docker` needs to be installed and the permissions for `docker`
must be given. Then, the deployment process should be automated like this:

```
./bin/deploy.sh
```
