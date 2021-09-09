import hammingdist
import numpy as np


def write_fasta_file(filename):
    with open(filename, "w") as f:
        f.write(">seq0\n")
        f.write("ACGTGTCGTGTCGACGTGTCG\n")
        f.write(">seq1\n")
        f.write("ACGTGTCGTTTCGACGAGTCG\n")
        f.write(">seq2\n")
        f.write("ACGTGTCGTTTCGACGAGTCG\n")
        f.write(">seq3\n")
        f.write("ACGTGTCGTTTCGACGAGTCG\n")
        f.write(">seq4\n")
        f.write("ACGTGTCGTGTCGACGTGTCG\n")
        f.write(">seq5\n")
        f.write("ACGTGTCGTATCGACGTGTCG\n")


def check_output_size(dat, n_in, n_out):
    tmp_out_file = "out.txt"

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


def test_from_fasta():
    tmp_in_file = "in.txt"
    write_fasta_file(tmp_in_file)

    data = hammingdist.from_fasta(tmp_in_file)
    check_output_size(data, 6, 6)

    data = hammingdist.from_fasta(tmp_in_file, n=5)
    check_output_size(data, 5, 5)

    data = hammingdist.from_fasta(tmp_in_file, include_x=True)
    check_output_size(data, 6, 6)

    data = hammingdist.from_fasta(tmp_in_file, remove_duplicates=True)
    check_output_size(data, 6, 3)

    data = hammingdist.from_fasta(tmp_in_file, include_x=True, n=2)
    check_output_size(data, 2, 2)

    data = hammingdist.from_fasta(tmp_in_file, remove_duplicates=True, n=3)
    check_output_size(data, 3, 2)

    data = hammingdist.from_fasta(
        tmp_in_file, remove_duplicates=True, n=5, include_x=True
    )
    check_output_size(data, 5, 2)
