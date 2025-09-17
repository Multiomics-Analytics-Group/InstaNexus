import os
import sys
import unittest

import scaffolding

script_dir = os.getcwd()  # get the current working directory
sys.path.append(os.path.join(script_dir, "../src"))


class TestOverlap(unittest.TestCase):

    def test_no_overlap(self):
        list1 = ["AAAA", "BBBB", "CCCC", "DDDD"]
        self.assertEqual(scaffolding.find_overlaps(list1, 2), [])

    def test_overlap_two_contigs(self):
        list1 = ["ABCD", "CDEF"]
        self.assertEqual(scaffolding.find_overlaps(list1, 2), [("ABCD", "CDEF", 2)])

    def test_overlap_three_contigs(self):
        list1 = ["ABCD", "CDEF", "EFGH"]
        self.assertEqual(
            scaffolding.find_overlaps(list1, 2),
            [("ABCD", "CDEF", 2), ("CDEF", "EFGH", 2)],
        )


class TestFilterContainedSequences(unittest.TestCase):

    def no_sequence_contained(self):
        list1 = ["ABCDEFGHI", "LMNO"]
        self.assertEqual(
            scaffolding.filter_contained_sequences(list1), ["ABCDEFGHI", "LMNO"]
        )


if __name__ == "__main__":
    # unittest.main() # this will run all the tests

    # unittest.main(defaultTest='TestOverlap.test_no_overlap')
    unittest.main(defaultTest="TestFilterContainedSequences.no_sequence_contained")
