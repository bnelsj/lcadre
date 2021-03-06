#! /usr/bin/env python3

"""
LCaDRE: Library Complexity and Duplication Rate Estimation
"""

__author__ = "Brad Nelson"
__contact__ = "bnelsj@gmail.com"

import argparse
from collections import Counter
import logging
from math import exp
import os

import pysam

ALLOWED_FILE_TYPES = set(["sam", "bam", "cram"])
MODE_DICT = {
    "sam": "r",
    "bam": "rb",
    "cram": "rc"
}

format = "[%(name)s - %(asctime)s] %(message)s"
logging.basicConfig(format=format, level=logging.INFO)
logger = logging.getLogger("lcadre")


class Lcadre(object):
    """
    Class for estimating library complexity and duplication rate from an alignment file.
    """
    def __init__(self, alignment_file, n_extrapolation, file_type="infer", **kwargs):
        self.alignment_file = alignment_file
        self.n_extrapolation = n_extrapolation
        self.file_type = file_type

        if self.file_type == "infer":
            ft = os.path.splitext(self.alignment_file)[1].lstrip(".")
            if ft not in ALLOWED_FILE_TYPES:
                msg = f"File type must be in {ALLOWED_FILE_TYPES}, instead it is {self.file_type}"
                logger.error(msg)
                raise ValueError(msg)
            else:
                self.file_type = ft
        elif self.file_type not in ALLOWED_FILE_TYPES:
            msg = f"File type must be in {ALLOWED_FILE_TYPES}, instead it is {self.file_type}"
            logger.error(msg)
            raise ValueError(msg)
        self.open_mode = MODE_DICT[self.file_type]
        self.counter = Counter()
        self.r_cutoff = 10 # Recommendation from Chao and Shen 2004

    def calc_s_ind(self, f0_hat):
        """Estimate number of signatures at target read pairs. Implements equation 9 from
        https://doi.org/10.1093/jpe/rtr044
        Arguments:
            f0_hat (float): Estimated number of undetected signatures
        Uses:
            self.m_star (int): Number of additional read pairs to collect
            self.total_read_pairs (int): Number of total read pairs in sample
            self.singletons (int): Count of signatures observed exactly once in sample
            self.total_signatures (int): Total distinct read signatures found
        Returns:
            s_ind (float): Estimated number of distinct signatures at target read pair count
        """
        expression = 1 - exp(-1 * self.m_star / self.total_read_pairs * self.singletons / f0_hat)
        s_ind = self.total_signatures + f0_hat * expression
        return s_ind

    def get_pos(self, read):
        """Get 5' start coordinate of read and correct for soft clipping.
        Arguments:
            read (pysam.AlignedSegment): Sequence read to get position for
        Returns:
            pos (int or None): 5' start coordinate of read, or None if unmapped
        """
        if read.reference_start is None or read.reference_end is None:
            pos = None
        elif read.is_reverse:
            pos = read.reference_end + read.query_alignment_start
        else:
            pos = read.reference_start - read.query_alignment_start
        return pos

    def get_read_string(self, read, pos):
        """Create a key for a read by finding the 5' start of the read, accounting for soft clipping.
        Arguments:
            read (pysam.AlignedSegment): AlignedSegment to get read string for
            pos (int or None): 5' start of read
        Returns:
            read_string (str): string to be used for key in counter
        """
        if read.reference_start is None or read.reference_end is None:
            return "-1"
        elif read.is_reverse:
            strand = "-"
        else:
            strand = "+"
        return f"{read.reference_name}.{pos}.{strand}"

    def order_pair(self, read_a, pos_a, read_b, pos_b):
        """Order reads in pair based on mapping, reference, corrected position, and strand.
        Arguments:
            read_a (pysam.AlignedSegment): First read in pair
            pos_a (int or None): 5' position of first read
            read_b (pysam.AlignedSegment): Second read in pair
            pos_b (int or None): 5' position of second read
        Returns:
            read_a, pos_a, read_b, pos_b: Reordered input based on mapping position (a and b may be flipped)
        """
        flip = False
        if read_a.reference_name is None and read_b.reference_name is None:
            pass
        elif read_b.reference_name is None:
            pass
        elif read_a.reference_name is None and read_b.reference_name is not None:
            flip = True
        elif read_a.reference_name > read_b.reference_name:
            flip = True
        elif read_a.reference_name == read_b.reference_name:
            if pos_a is None and pos_b is not None:
                flip = True
            elif pos_a is not None and pos_b is not None and pos_a > pos_b:
                flip = True
        if flip:
            read_b, read_a = read_a, read_b
            pos_b, pos_a = pos_a, pos_b

        return read_a, pos_a, read_b, pos_b

    def pair_generator(self):
        """Pair alignment records and yield them with a generator.

        Uses:
            self.alignment_file (str): Path to alignment file
            self.open_mode (str): Open mode for alignment file (See MODE_DICT)
        Yields:
            (read_a, read_b) (pysam.AlignedSegment, pysam.AlignedSegment): Primary alignments for next read pair
        """
        with pysam.AlignmentFile(self.alignment_file, self.open_mode) as reader:
            read_a = None
            read_b = None

            for read in reader:
                if read.is_secondary or read.is_supplementary:
                    continue
                elif read_a is None:
                    read_a = read
                elif read.query_name == read_a.query_name:
                    read_b = read
                    yield (read_a, read_b)
                    read_a = None
                    read_b = None
                else:
                    read_a = read

    def perform_estimation(self):
        """
        Estimate library complexity using the Chao1  and ACE diversity estimators
        and extrapolate the duplicate rate using sampling theory.
        See https://doi.org/10.1093/jpe/rtr044 for a review.

        Uses:
            self.counter (collections.Counter): Modified dict with counts for each signature
            self.r_cutoff (int): Cutoff between rare and common signatures
        Sets:
            self.total_read_pairs (int): Total read pairs in primary alignments
            self.total_signatures (int): Total distinct read signatures found
            self.singletons (int): Count of signatures found exactly once
            self.doubletons (int): Count of signatures found exactly twice
            self.chao1 (float): Lower estimate of total signatures present in library
            self.m_star (int): Number of additional reads to collect
            self.f0_hat (float): Estimate of undiscovered signatures remaining in library
            self.s_ind (float): Estimate of signatures with self.n_extrapolation read pairs sequenced
            self.dup_rate_obs (float): Observed duplicate rate
            self.dup_rate_extrap (float): Extrapolated duplicate rate at self.n_extrapolation read pairs

            self.s_rare (int): Count of rare signatures
            self.x_rare (int): Abundance of rare signatures
            self.k_counter (Counter): Frequency count for values less than r_cutoff
            self.c_ace (float): Sample coverage estimate
            self.cov_var_sq (float): Sample coverage variance (squared)
            self.ace (float): Estimate of undiscovered signatures remaining in library (ACE)
            self.s_ind_ace (float): Estimate of signatures at target read pair count (ACE)
        """
        self.singletons = 0
        self.doubletons = 0
        others = 0

        self.s_rare = 0 # Count of rare signatures
        self.x_rare = 0 # Abundance of rare signatures
        self.k_counter = Counter() # Count of signatures for each freq at or below r cutoff

        for pos, count in self.counter.items():
            if count <= self.r_cutoff:
                self.s_rare += 1
                self.x_rare += count
                self.k_counter[count] += 1
            if count == 1:
                self.singletons += 1
            elif count == 2:
                self.doubletons += 1
            else:
                others += 1

        if self.doubletons == 0:
            msg = "No doubletons found. Likely due to too few reads."
            logger.error(msg)
            raise ZeroDivisionError(msg)

        self.total_read_pairs = sum(self.counter.values())
        self.total_signatures = self.singletons + self.doubletons + others

        self.chao1 = self.total_signatures + self.singletons ** 2 / self.doubletons / 2

        # Number of additional reads to collect
        self.m_star = self.n_extrapolation - self.total_read_pairs

        # Estimate of remaining unique signatures in library
        self.f0_hat = self.chao1 - self.total_signatures

        self.s_ind = self.calc_s_ind(self.f0_hat)

        self.dup_rate_obs = 1 - self.total_signatures / self.total_read_pairs
        self.dup_rate_extrap = 1 - self.s_ind / self.n_extrapolation

        # Use Abundance-based Coverage Estimator (ACE) to improve estimation of remaining total_signatures
        # Equation 16 from https://doi.org/10.1093/jpe/rtr044

        self.c_ace = 1 - self.singletons / self.x_rare # Sample coverage estimate

        top_sum = sum(k * (k - 1) * self.k_counter[k] for k in range(1, self.r_cutoff + 1))
        bot_sum = sum(k * self.k_counter[k] * (k * self.k_counter[k] - 1) for k in range(1, self.r_cutoff + 1))
        self.cov_var_sq = max(self.s_rare / self.c_ace * top_sum / bot_sum - 1, 0)
        self.f0_hat_ace = self.s_rare / self.c_ace + self.singletons / self.c_ace * self.cov_var_sq - self.s_rare
        self.ace = self.f0_hat_ace + self.total_signatures

        self.s_ind_ace = self.calc_s_ind(self.f0_hat_ace)
        self.dup_rate_extrap_ace = 1 - self.s_ind_ace / self.n_extrapolation

    def process_pair(self, read_pair):
        """Process a read pair and add it to the counter.
        Arguments:
            read_pair (2-tuple(pysam.AlignedSegment)): Pair of reads that are primary alignments
        Uses:
            self.counter (collections.Counter): Counter for number of times read pair signature has been seen.
        """
        read_a, read_b = read_pair
        pos_a = self.get_pos(read_a)
        pos_b = self.get_pos(read_b)
        read_a, pos_a, read_b, pos_b = self.order_pair(read_a, pos_a, read_b, pos_b)

        str_a = self.get_read_string(read_a, pos_a)
        str_b = self.get_read_string(read_b, pos_b)

        # Key signature inspired by https://github.com/GregoryFaust/samblaster#DupIdentification
        key = f"{str_a}.{str_b}"

        self.counter[key] += 1

    def read_input(self):
        """Read alignment file and process all found read pairs
        """
        logger.info(f"Reading alignment file {self.alignment_file}")
        pair_generator = self.pair_generator()
        for read_pair in pair_generator:
            self.process_pair(read_pair)

    def report(self):
        """Log output of estimation procedure.
        Uses:
            self.total_read_pairs (int): Total read pairs found
            self.total_signatures (int): Total distinct read signatures found
            self.singletons (int): Count of signatures found exactly once
            self.doubletons (int): Count of signatures found exactly twice
            self.dup_rate_obs (float): Observed duplicate rate
            self.n_extrapolation (int): Total read pairs target for extrapolation
            self.chao1 (float): Lower estimate of total signatures present in library
            self.s_ind (float): Estimated distinct signatures at read pair target using Chao1
            self.dup_rate_extrap (float): Estimated duplicate rate at target using Chao1
            self.s_ind_ace (float): Estimated distinct signatures at read pair target using ACE
            self.ace (float): Estimate of distinct signatures present in library using ACE
            self.dup_rate_extrap_ace (float): Estimated duplicate rate at target using ACE

        """
        logger.info(f"Found {self.total_read_pairs:,} read pairs")
        logger.info(f"Found {self.total_signatures:,} non-duplicate signatures")
        logger.info(f"Found {self.singletons:,} signatures that occurred exactly once")
        logger.info(f"Found {self.doubletons:,} signatures that occurred exactly twice")
        logger.info(f"Observed duplicate rate: {self.dup_rate_obs:0.3f}")
        logger.info("")
        logger.info(f"##### Chao1-based estimates #####")
        logger.info(f"Library complexity estimate: {self.chao1:,.0f}")
        logger.info(f"Estimated signatures at {self.n_extrapolation:,} read pairs: {self.s_ind:,.0f}")
        logger.info(f"Extrapolated duplication rate: {self.dup_rate_extrap:0.3f}")
        logger.info("")
        logger.info(f"##### ACE-based estimates #####")
        logger.info(f"Estimated signatures at {self.n_extrapolation:,} read pairs: {self.s_ind_ace:,.0f}")
        logger.info(f"Library complexity estimate: {self.ace:,.0f}")
        logger.info(f"Extrapolated duplication rate: {self.dup_rate_extrap_ace:0.3f}")

    def run(self):
        """Run LCaDRE analysis.
        """
        self.verify_header()
        self.read_input()
        self.perform_estimation()
        self.report()

    def verify_header(self):
        """Confirm the alignment file is compatible with LCaDRE using header metadata.
        Uses:
            self.alignment_file (str): Path to alignment file
            self.open_mode (str): Open mode for alignment file (See MODE_DICT)
        Raises:
            ValueError: if header indicates file is coordinate sorted
        """
        with pysam.AlignmentFile(self.alignment_file, self.open_mode) as reader:
            header = reader.header

        if "HD" in header and "SO" in header["HD"] and header["HD"]["SO"].lower().strip() == "coordinate":
            msg = "lcadre.py requires read name sorted alignment files, but the file is coordinate sorted"
            logger.error(msg)
            raise ValueError(msg)


def parse_args():
    parser = argparse.ArgumentParser(description="LCaDRE: Library Complexity and Duplication Rate Estimation")
    parser.add_argument("alignment_file", help="Path to read name-sorted alignment file")
    parser.add_argument("--file_type",
                        "-t",
                        default="infer",
                        choices=["sam", "bam", "cram", "infer"],
                        help="Type of alignment file (Default: %(default)s)")
    parser.add_argument("--n_extrapolation", "-n", default=100000000, type=int,
                        help="Target count of read pairs for extrapolation (Default: %(default)s)")

    args = vars(parser.parse_args())

    return args


def main():
    args = parse_args()
    logger.info("LCaDRE: Library Complexity and Duplication Rate Estimation")
    lcadre = Lcadre(**args)
    lcadre.run()


if __name__ == "__main__":
    main()
