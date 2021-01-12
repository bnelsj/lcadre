[![Build Status](https://travis-ci.com/bnelsj/lcadre.svg?branch=main)](https://travis-ci.com/bnelsj/lcadre)

# LCaDRE: Library Complexity and Duplication Rate Estimation

## Description

When we want to know if a new DNA library preparation was successful, it is often useful to sequence it a little and measure performance. This can be done for many metrics with around 500K read pairs, but for duplicate rate, it is more difficult.

For a particular sequencing run, duplicate rate could be driven by a few very common read pair signatures, in which case duplicate rate may not change significantly in a full run, or it could be driven by low library complexity, in which case deeper sequencing can lead to saturation and much higher duplicate rates.

This project aims to estimate library complexity and extrapolate duplication rate from the given sequence data to a full sequencing run using methods from sampling theory. LCaDRE uses a read pair signature definition inspired by samblaster, but instead of just flagging duplicates, it records how many times each signature has been seen. It can then use the Chao1 diversity estimator (a non-parametric lower bound estimator for species diversity) to estimate the library complexity. To extrapolate the duplication rate to a larger sequencing run, it uses equations from [Colwell et al. 2012](https://doi.org/10.1093/jpe/rtr044). While these methods were developed for species diversity estimation in ecology, much of the early work traces back to [Alan Turing's work cracking the Enigma machine](https://en.wikipedia.org/wiki/Good%E2%80%93Turing_frequency_estimation).

## Installation

Create the conda environment using the included env.yml file:

```
conda env create -f env.yml
conda activate lcadre
```

## Usage

```
python lcadre.py <alignment file> -n <read pairs in full run>
```

## Tests

Tests can be run with the following command:
```
python -m unittest discover
```

## Example output

```
[lcadre - 2021-01-11 17:49:49,279] LCaDRE: Library Complexity and Duplication Rate Estimation
[lcadre - 2021-01-11 17:49:49,280] Reading alignment file NA12878.rn.bam
[lcadre - 2021-01-11 17:50:13,705] Found 3,685,262 read pairs
[lcadre - 2021-01-11 17:50:13,705] Found 3,508,436 distinct signatures
[lcadre - 2021-01-11 17:50:13,705] Found 3,338,134 signatures that occurred exactly once
[lcadre - 2021-01-11 17:50:13,705] Found 163,999 signatures that occurred exactly twice
[lcadre - 2021-01-11 17:50:13,705] Estimated distinct signatures remaining in library: 33,973,191
[lcadre - 2021-01-11 17:50:13,705] Library complexity estimate: 37,481,627
[lcadre - 2021-01-11 17:50:13,705] Estimated signatures at 100,000,000 read pairs: 34,876,215
[lcadre - 2021-01-11 17:50:13,705] Observed duplication rate: 0.048
[lcadre - 2021-01-11 17:50:13,706] Extrapolated duplication rate: 0.651
```

## Evaluation

Tests were performed on chr11 of the low coverage 1000 genomes dataset for NA12878, available [here](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12878/alignment/NA12878.chrom11.ILLUMINA.bwa.CEU.low_coverage.20121211.bam). Picard's EstimateLibraryComplexity tool was used for comparison.

| | LCaDRE | Picard ELC |
| --- | --- | --- |
| Observed duplication rate | 0.048 | 0.056 |
| Estimated library complexity | 37481627 | 27065571 |
| Extrapolated duplication rate (100M read pairs) | 0.651 | NA |
| Runtime (seconds) | 25 | 31 |

It's good that the observed duplication rate is so close to Picard's measure, and it's not surprising that the library complexity estimates are so different--this kind of extrapolation is imprecise, and confidence intervals would likely overlap. I was initially surprised to see that lcadre's runtime is faster than Picard's, since python tends to be much slower than java, but picard is doing a lot of additional things that can be somewhat time consuming, including actually looking at query sequence.

At this point, I would trust picard's library complexity estimate over lcadre's due to the former's maturity and wide use. I was unable to find public libraries with high duplicate rate, which would be extremely useful for testing the accuracy of the extrapolation under more conditions.
