import hammingdist
import numpy as np
import random
import pytest

gpu_options = [False, True] if hammingdist.cuda_gpu_available() else [False]


def write_fasta_file(filename, sequences):
    with open(filename, "w") as f:
        for i, seq in enumerate(sequences):
            f.write(f">seq{i}\n{seq}\n")


def check_output_sizes(dat, n_in, n_out, tmp_out_file, fasta_sequence_indices=None):
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
    if fasta_sequence_indices is not None:
        assert np.allclose(indices, fasta_sequence_indices)


@pytest.mark.parametrize(
    "from_fasta_func", [hammingdist.from_fasta, hammingdist.from_fasta_large]
)
@pytest.mark.parametrize("use_gpu", gpu_options)
def test_from_fasta(from_fasta_func, use_gpu, tmp_path):
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

    data = from_fasta_func(fasta_file, use_gpu=use_gpu)
    check_output_sizes(data, 6, 6, output_file)

    data = from_fasta_func(fasta_file, n=5, use_gpu=use_gpu)
    check_output_sizes(data, 5, 5, output_file)

    data = from_fasta_func(fasta_file, include_x=True, use_gpu=False)
    check_output_sizes(data, 6, 6, output_file)

    fasta_sequence_indices = hammingdist.fasta_sequence_indices(fasta_file)
    data = from_fasta_func(fasta_file, remove_duplicates=True, use_gpu=use_gpu)
    check_output_sizes(data, 6, 4, output_file, fasta_sequence_indices)

    fasta_sequence_indices = hammingdist.fasta_sequence_indices(fasta_file)
    data = from_fasta_func(
        fasta_file, remove_duplicates=True, include_x=True, use_gpu=False
    )
    check_output_sizes(data, 6, 4, output_file, fasta_sequence_indices)

    fasta_sequence_indices = hammingdist.fasta_sequence_indices(fasta_file, n=2)
    data = from_fasta_func(fasta_file, include_x=True, n=2, use_gpu=False)
    check_output_sizes(data, 2, 2, output_file, fasta_sequence_indices)

    fasta_sequence_indices = hammingdist.fasta_sequence_indices(fasta_file, n=3)
    data = from_fasta_func(fasta_file, remove_duplicates=True, n=3, use_gpu=use_gpu)
    check_output_sizes(data, 3, 2, output_file, fasta_sequence_indices)

    fasta_sequence_indices = hammingdist.fasta_sequence_indices(fasta_file, n=5)
    data = from_fasta_func(
        fasta_file, remove_duplicates=True, n=5, include_x=True, use_gpu=False
    )
    check_output_sizes(data, 5, 3, output_file, fasta_sequence_indices)


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
@pytest.mark.parametrize("max_distance", [0, 1, 2, 3, 89, 497, 9999999])
def test_fasta_reference_distances(chars, include_x, max_distance, tmp_path):
    # generate 50 sequences, each with 25 characters
    sequences = ["".join(random.choices(chars, k=25)) for i in range(50)]
    fasta_file = str(tmp_path / "fasta.txt")
    write_fasta_file(fasta_file, sequences)
    # calculate distances matrix
    data = hammingdist.from_fasta(
        fasta_file,
        remove_duplicates=False,
        include_x=include_x,
        max_distance=max_distance,
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
                assert data[i, j] == min(max_distance, dist)
            # should also agree with output of distance function for these two sequences
            assert dist == hammingdist.distance(
                sequences[i], sequences[j], include_x=include_x
            )


def test_distance():
    assert hammingdist.distance("ACGT", "ACCT") == 1
    # here X is invalid so has distance 1 from itself:
    assert hammingdist.distance("ACGTX", "ACCTX") == 2
    # now X is valid so has distance 0 from itself:
    assert hammingdist.distance("ACGTX", "ACCTX", include_x=True) == 1
    with pytest.raises(RuntimeError):
        hammingdist.distance("ACGT", "ACC")


@pytest.mark.skipif(
    not hammingdist.cuda_gpu_available(),
    reason="No CUDA GPU available or hammingdist was compiled without CUDA support",
)
@pytest.mark.parametrize("max_distance", [0, 1, 2, 3, 89, 497, 9999999])
def test_from_fasta_to_lower_triangular(tmp_path, max_distance):
    sequences = [
        "ACGTGTCGTGTCGACGTGTCGCAGGTGTCGACGTGTCGCAGGTGTCGACGTGTCGCAG",
        "CCGTGTCGTGTCGACGTGTCGC-GGTGTCGACGTGTCGCAGGTGTCGACGTGTCGCAG",
        "CAGTGT-GTGTCGACGTGTCGCAGGTGTCGACGTGTCGCAGGTGTCGACGTG--GCAG",
    ]
    lower_triangular_dist = [
        [min(1, max_distance)],
        [min(2, max_distance), min(1, max_distance)],
    ]
    fasta_file = str(tmp_path / "fasta.txt")
    output_file = str(tmp_path / "out.txt")
    write_fasta_file(fasta_file, sequences)
    hammingdist.from_fasta_to_lower_triangular(
        fasta_file, output_file, max_distance=max_distance
    )
    with open(output_file) as f:
        data = f.read().splitlines()
    assert len(data) == 2
    assert np.allclose(np.fromstring(data[0], sep=","), lower_triangular_dist[0])
    assert np.allclose(np.fromstring(data[1], sep=","), lower_triangular_dist[1])


@pytest.mark.parametrize("samples", [2, 3, 5, 11, 54, 120, 532, 981, 1568])
@pytest.mark.parametrize("threshold", [0, 1, 2, 3, 4, 9, 89, 497])
def test_dump_sparse(tmp_path, samples, threshold):
    lt_file = str(tmp_path / "lt.txt")
    with open(lt_file, "w") as f:
        for n in range(1, samples):
            f.write(
                ",".join(
                    str(x)
                    for x in np.random.randint(low=0, high=256, size=n, dtype=np.uint32)
                )
            )
            f.write("\n")
    data = hammingdist.from_lower_triangular(lt_file)
    # expected number of sparse entries
    num_below_threshold = np.sum(np.array(data._distances) <= threshold)
    if num_below_threshold > 0:
        sparse_file = str(tmp_path / "sparse.txt")
        data.dump_sparse(sparse_file, threshold)
        sparse = np.loadtxt(sparse_file, delimiter=" ", dtype=np.uint32)
        if sparse.ndim == 1:
            sparse = np.array([sparse])
        assert sparse.shape[0] == num_below_threshold
        assert sparse.shape[1] == 3
        # each sparse entry has the correct distance, and is not above threshold
        for e in sparse:
            assert e[2] <= threshold
            assert data[e[0], e[1]] == e[2]
