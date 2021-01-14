#!usr/bin/env/python3

from collections import Counter
import logging
import os
import unittest

import lcadre

logger = logging.getLogger("lcadre")
logger.setLevel(logging.ERROR)


class MockPair:
    def __init__(self,
                 reference_name=None,
                 reference_start=None,
                 reference_end=None,
                 is_reverse=False,
                 query_alignment_start=0):
        self.reference_name = reference_name
        self.reference_start = reference_start
        self.reference_end = reference_end
        self.is_reverse = is_reverse
        self.query_alignment_start = query_alignment_start

    def __eq__(self, other):
        return self.reference_name == other.reference_name and \
               self.reference_start == other.reference_start and \
               self.reference_end == other.reference_end and \
               (self.is_reverse == other.is_reverse) and \
               self.query_alignment_start == other.query_alignment_start

    def __repr__(self):
        return f"{self.reference_name}:{self.reference_start}-{self.reference_end}:{self.is_reverse}:{self.query_alignment_start}"

class TestLcadre(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        pass

    def setUp(self):

        self.dirname = os.path.dirname(__file__)
        self.collateral_dir = self.dirname + "/collateral/"
        self.input_dir = self.collateral_dir + "input/"

        self.n = 100000000
        self.fn = os.path.join(self.input_dir, "simple.bam")
        self.empty_fn = os.path.join(self.input_dir, "empty.bam")

        self.test = lcadre.Lcadre(self.fn, self.n)

        self.read_blank = MockPair()
        self.read_a = MockPair(reference_name="chr1", reference_start=100, reference_end=200, query_alignment_start=10)
        self.read_b = MockPair(reference_name="chr1", reference_start=100, reference_end=200, query_alignment_start=10, is_reverse=True)
        self.read_c = MockPair(reference_name="chr2", reference_start=100, reference_end=200, query_alignment_start=10)

    def tearDown(self):
        pass

    def testFileTypes(self):
        fts = ["sam", "bam", "cram", "infer"]

        for ft in fts:
            test = lcadre.Lcadre(self.fn, self.n, ft)
            self.assertTrue(test.file_type in lcadre.ALLOWED_FILE_TYPES)

        fail_types = ["fasta", "", None, 1]

        for ft in fail_types:
            with self.assertRaises(ValueError):
                test = lcadre.Lcadre(self.fn, self.n, ft)

        with self.assertRaises(ValueError):
            test = lcadre.Lcadre("test.txt", self.n)

    def testEstimation(self):
        counter = Counter()

        for item in [1, 1, 2, 4, 1, 3, 4, 5, 0]:
            counter[item] += 1

        test = lcadre.Lcadre(self.fn, self.n)
        test.counter = counter
        test.perform_estimation()

        self.assertEqual(test.singletons, 4)
        self.assertEqual(test.doubletons, 1)
        self.assertEqual(test.total_read_pairs, 9)
        self.assertEqual(test.total_signatures, 6)
        self.assertEqual(test.chao1, 14)
        self.assertEqual(test.m_star, self.n - 9)
        self.assertEqual(test.f0_hat, 8)
        self.assertEqual(test.s_ind, 14)
        self.assertEqual(test.dup_rate_obs, 0.33333333333333337)
        self.assertEqual(test.dup_rate_extrap, 0.99999986)

    def testProcessPair(self):
        test = lcadre.Lcadre(self.fn, self.n)
        test.process_pair((self.read_blank, self.read_blank))
        self.assertTrue("-1.-1" in test.counter)
        self.assertEqual(test.counter["-1.-1"], 1)
        self.assertEqual(sum(test.counter.values()), 1)

    def testGetPos(self):
        self.assertEqual(self.test.get_pos(self.read_blank), None)
        self.assertEqual(self.test.get_pos(self.read_a), 90)
        self.assertEqual(self.test.get_pos(self.read_b), 210)

    def testGetReadString(self):
        self.assertEqual(self.test.get_read_string(self.read_blank, -1), "-1")
        self.assertEqual(self.test.get_read_string(self.read_a, 90), "chr1.90.+")
        self.assertEqual(self.test.get_read_string(self.read_b, 210), "chr1.210.-")

    def testOrderPair(self):
        read_a, pos_a, read_b, pos_b = self.test.order_pair(self.read_a, 90, self.read_b, 210)
        self.assertEqual(self.read_a, read_a)
        self.assertEqual(pos_a, 90)
        self.assertEqual(self.read_b, read_b)
        self.assertEqual(pos_b, 210)

        read_a, pos_a, read_b, pos_b = self.test.order_pair(self.read_b, 210, self.read_a, 90)
        self.assertEqual(self.read_a, read_a)
        self.assertEqual(pos_a, 90)
        self.assertEqual(self.read_b, read_b)
        self.assertEqual(pos_b, 210)

        read_b, pos_b, read_blank, pos_blank = self.test.order_pair(self.read_b, 210, self.read_blank, -1)
        self.assertEqual(self.read_b, read_b)
        self.assertEqual(pos_b, 210)
        self.assertEqual(self.read_blank, read_blank)
        self.assertEqual(pos_blank, -1)

        read_b, pos_b, read_blank, pos_blank = self.test.order_pair(self.read_blank, -1, self.read_b, 210)
        self.assertEqual(self.read_b, read_b)
        self.assertEqual(pos_b, 210)
        self.assertEqual(self.read_blank, read_blank)
        self.assertEqual(pos_blank, -1)

        read_b, pos_b, read_c, pos_c = self.test.order_pair(self.read_c, 90, self.read_b, 210)
        self.assertEqual(self.read_b, read_b)
        self.assertEqual(pos_b, 210)
        self.assertEqual(self.read_c, read_c)
        self.assertEqual(pos_c, 90)

        read_b, pos_b, read_a, pos_a = self.test.order_pair(self.read_a, None, self.read_b, 210)
        self.assertEqual(self.read_a, read_a)
        self.assertEqual(pos_a, None)
        self.assertEqual(self.read_b, read_b)
        self.assertEqual(pos_b, 210)

    def testRun(self):
        self.test.run()

    def testEmptyBam(self):
        test = lcadre.Lcadre(self.empty_fn, self.n)
        with self.assertRaises(ZeroDivisionError):
            test.run()


if __name__ == '__main__':
    unittest.main()
