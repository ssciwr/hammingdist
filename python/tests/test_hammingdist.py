import hammingdist
import numpy as np
import random
import pytest


def write_fasta_file(filename, sequences):
    with open(filename, "w") as f:
        for i, seq in enumerate(sequences):
            f.write(f">seq{i}\n{seq}\n")


def check_output_sizes(dat, n_in, n_out, tmp_out_file):
    dat.dump(tmp_out_file)
    dump = np.loadtxt(tmp_out_file, delimiter=",")
    assert len(dump) == n_out
    assert len(dump[0]) == n_out

    dat.dump_lower_triangular(tmp_out_file)
    with open(tmp_out_file) as f:
        data = f.read().splitlines()
    lt_first = np.fromstring(data[0], sep=",")
    lt_last = np.fromstring(data[-1], sep=",")
    assert len(lt_first) == 1
    assert len(lt_last) == n_out - 1

    dat.dump_sequence_indices(tmp_out_file)
    indices = np.loadtxt(tmp_out_file, delimiter=",")
    assert len(indices) == n_in
    assert indices[0] == 0


def test_from_fasta(tmp_path):
    sequences = [
        "ACGTGTCGTGTCGACGTGTCG",
        "ACGTGTCGTTTCGACGAGTCG",
        "ACGTGTCGTTTCGACGAGTCG",
        "ACGTGTCGTTTCGACGAGTCG",
        "ACGTGACGTGTCGACGTGTCG",
        "ACGTGXCGTGTCGACGTGTCG",
    ]
    fasta_file = str(tmp_path / "fasta.txt")
    output_file = str(tmp_path / "out.txt")
    write_fasta_file(fasta_file, sequences)

    data = hammingdist.from_fasta(fasta_file)
    check_output_sizes(data, 6, 6, output_file)

    data = hammingdist.from_fasta(fasta_file, n=5)
    check_output_sizes(data, 5, 5, output_file)

    data = hammingdist.from_fasta(fasta_file, include_x=True)
    check_output_sizes(data, 6, 6, output_file)

    data = hammingdist.from_fasta(fasta_file, remove_duplicates=True)
    check_output_sizes(data, 6, 4, output_file)

    data = hammingdist.from_fasta(fasta_file, remove_duplicates=True, include_x=True)
    check_output_sizes(data, 6, 4, output_file)

    data = hammingdist.from_fasta(fasta_file, include_x=True, n=2)
    check_output_sizes(data, 2, 2, output_file)

    data = hammingdist.from_fasta(fasta_file, remove_duplicates=True, n=3)
    check_output_sizes(data, 3, 2, output_file)

    data = hammingdist.from_fasta(
        fasta_file, remove_duplicates=True, n=5, include_x=True
    )
    check_output_sizes(data, 5, 3, output_file)


@pytest.mark.parametrize(
    "chars,include_x",
    [
        (["A"], False),
        (["A", "C", "G", "T", "-"], False),
        (["A", "C", "G", "T", "-"], True),
        (["A", "C", "G", "T", "-", "X"], False),
        (["A", "C", "G", "T", "-", "X"], True),
    ],
)
def test_fasta_reference_distances(chars, include_x, tmp_path):
    # generate 50 sequences, each with 25 characters
    sequences = ["".join(random.choices(chars, k=25)) for i in range(50)]
    fasta_file = str(tmp_path / "fasta.txt")
    write_fasta_file(fasta_file, sequences)
    # calculate distances matrix
    data = hammingdist.from_fasta(
        fasta_file, remove_duplicates=False, include_x=include_x
    )
    # use each sequence in turn as the reference sequence & calculate reference distances
    for i, sequence in enumerate(sequences):
        vec = hammingdist.fasta_reference_distances(
            sequence, fasta_file, include_x=include_x
        )
        assert len(vec) == len(sequences)
        assert vec.dtype == np.uint32
        # reference distance vector should match corresponding column of distances matrix
        for j, dist in enumerate(vec):
            # if x is not included, invalid chars have distance 1 but data[i,i] returns 0 by construction
            if include_x or i != j:
                assert data[i, j] == dist
