import unittest

import greedy_method


class TestOverlap(
    unittest.TestCase
):  # testcase class indicates that this is a test case

    def test_no_overlap(self):

        peptides1 = ["AAAA", "BBBB", "CCCC", "DDDD"]
        self.assertEqual(greedy_method.find_peptide_overlaps(peptides1, 2), {})

    def test_single_overlap(self):

        peptides1 = ["ABCD", "CDEF"]
        self.assertEqual(
            greedy_method.find_peptide_overlaps(peptides1, 2), {0: [(1, 2), (1, 2)]}
        )

        peptides2 = ["ABCD", "EFAB"]
        self.assertEqual(
            greedy_method.find_peptide_overlaps(peptides2, 2), {1: [(0, 2), (0, 2)]}
        )

        peptides3 = ["ABCDE", "CDEFG"]
        self.assertEqual(
            greedy_method.find_peptide_overlaps(peptides3, 2), {0: [(1, 3), (1, 3)]}
        )


class TestAssemble(unittest.TestCase):

    def test_two_peptides(self):

        peptides1 = ["ABCD", "CDEF"]
        self.assertEqual(
            greedy_method.assemble_contigs(peptides1, min_overlap=2), ["ABCDEF"]
        )

        peptides2 = ["ABCD", "EFAB"]
        self.assertEqual(
            greedy_method.assemble_contigs(peptides2, min_overlap=2), ["EFABCD"]
        )

        peptides3 = ["ABCDE", "CDEFG"]
        self.assertEqual(
            greedy_method.assemble_contigs(peptides3, min_overlap=2), ["ABCDEFG"]
        )

        peptides4 = ["ABCD", "CDEF", "EFGH"]
        self.assertEqual(
            greedy_method.assemble_contigs(peptides4, min_overlap=2), ["ABCDEFGH"]
        )

    def test_multiple_peptides(self):
        # exception case where the peptides are not overlapping
        peptides5 = ["ABCD", "CDEF", "ZZXX"]
        self.assertEqual(
            greedy_method.assemble_contigs(peptides5, min_overlap=2), ["ABCDEF", "ZZXX"]
        )

        # exception case where the peptides are not overlapping
        # SPECIAL CASE TO UNDERSTAND BETTER
        peptides6 = ["ABCD", "CDAB", "ZZXX"]
        self.assertEqual(
            greedy_method.assemble_contigs(peptides6, min_overlap=2), ["ABCDAB", "ZZXX"]
        )


class TestFindContigOverlap(unittest.TestCase):

    def test_overlap_found(self):
        # using two sequences with a clear overlap "CCC"
        # setting the minimum overlap to 2 even though the overlap is 3
        # since the sequences have equal length seq2 is treated as the smaller sequence.

        seq1 = "AAACCC"
        seq2 = "CCCGGG"
        min_overlap = 2

        result = greedy_method.find_contig_overlap(seq1, seq2, min_overlap)
        self.assertIsNotNone(result, "Expected to find a valid overlap, but got None.")
        # the string is printed if the result is None.

        (
            smallest_seq,
            largest_seq,
            max_overlap_len,
            max_overlap_pos_small,
            max_overlap_pos_large,
        ) = result

        self.assertEqual(largest_seq, seq1)  # seq1 is the larger sequence.
        self.assertEqual(
            smallest_seq, seq2
        )  # seq2 is the smaller sequence (in this case is considered smaller because the second one).
        self.assertEqual(
            max_overlap_len, 3
        )  # expected maximum overlap is "CCC" with length 3.
        self.assertEqual(
            max_overlap_pos_small, 0
        )  # seq2 ("CCCGGG") overlap "CCC" starts at position 0
        self.assertEqual(
            max_overlap_pos_large, 3
        )  # seq1 ("AAACCC") overlap "CCC" starts at position 3.

    def test_no_overlap(self):
        # test case where no valid overlap is expected.
        seq1 = "AAAAAA"
        seq2 = "TTTTTT"
        min_overlap = 2

        result = greedy_method.find_contig_overlap(seq1, seq2, min_overlap)
        self.assertIsNone(result, "Expected no overlap, but got a result.")


if __name__ == "__main__":
    # unittest.main() # this will run all the tests
    # unittest.main(defaultTest='TestOverlap.test_single_overlap')
    # unittest.main(defaultTest='TestAssemble.test_multiple_peptides')
    unittest.main(defaultTest="TestFindContigOverlap.test_no_overlap")
