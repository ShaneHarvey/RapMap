# Efficiently solve read mapping alignments

# RapMap TODOS
SAM flags field is wrong (0x900 for only first read mapping)?

Add alignment to paired-end reads. Should be easy!

Optimize RapMapAligner?

Group matricies (m, x, y, tm, tx, ty) for better cache preformance?

Cache transcript/read to avoid repeated alignments! This could lead to a big
performance increase because often the hit regions are exactly the same bewteen
transcripts.

Or could remember the Next Informative Position of the match (but what abut the
before the match part?) Also this is complicated because the QuasiAlignment hits
could be for different parts of the SA ie they are hits from a different k-mer.
This means that care would need to be taken so that we don't compare the NIP from
hits originating from different k-mers in the read. This still ignores the part of
the transcript that matches with the read before the k-mer.

QuasiAlignment hits would need to be grouped by the origniating SA interval hit
(could use the saved query position or *queryPos.*) Then we can see if the after
match alignment can be reused for each alignment in this group buy comparing the
saved NIP. This extra management should not be more complex than doing an extra
semi-global alignment.

Test cases to make sure aligner and CIGAR reporting works properly.


# Some alignments are sub-optimal!
Optimal Semi-global alignment:
Trans: TCCAGGGGGAGGCCTGCAGGCCCCTGGCCCCTTCCACCACCTCTGCCCTCCGTCTGCAGACCTCGTCCATCTGCACCAGGCTCTGCCTTCACTCCCCCAAGTCTTTGAAAATTTGTTCCTTTCCTTTGAAGTCACATTTTCTTTTAAAATTTTTTG
Read:  GATTCGCCTTGTGCCTTTATACCGACCTCATCTGCACTGGGCTCTGCCTTCACTCCCCCAAGTCTTTGAAAATTTA
Score: 18
Alignment matrix:
      0     .    :    .    :    .    :    .    :    .    :
        TCCAGGGGGAGGCCTGCAGGCCCCTGGCCCCTTCCACCACCTCTGCCCTC

        --------------------------------------------------

     50     .    :    .    :    .    :    .    :    .    :
        CG--TC------TGC----A----GACCTCGTCCATCTGCACCAGGCTCT
         |  ||      |||    |    |||||    |||||||||  ||||||
        -GATTCGCCTTGTGCCTTTATACCGACCT----CATCTGCACTGGGCTCT

    100     .    :    .    :    .    :    .    :    .    :
        GCCTTCACTCCCCCAAGTCTTTGAAAATTT-GTTCCTTTCCTTTGAAGTC
        ||||||||||||||||||||||||||||||
        GCCTTCACTCCCCCAAGTCTTTGAAAATTTA-------------------

    150     .    :    .    :
        ACATTTTCTTTTAAAATTTTTTG

        -----------------------

RapMap Sub-Optimal Alignment:
Score: 6
CIGAR: 1=4I2=2I1=1X6=1X1=1D1=1I2=2I2=1I2=2D2=2D2=1X1=3D1X2=2D36=1X



Trans: ACTTTTGTAAAGATTAAGCTCATTTAGTGT TGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT AGTATTTCAGCAGGATCTGCTGGCAGGGTTTTTTTGTTTTATTTGTTTGCTTATTTTTAAATTAACTGTTTTGAGCTTTGA
Read:  GCCTTCCGTGCCTTG                TGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT TTTTTTTTTTTTTTTTTTTTTT
Reverse hit!
Score: 11
Alignment matrix:
      0     .    :    .    :    .    :    .    :    .    :
        ACTTTTGTAAAGATTAAGCTCATTTAGTG--TTGT-TTTTTTTTTTTTTT
                         || | ||  |||  |||| ||||||||||||||
        -----------------GC-C-TTCCGTGCCTTGTGTTTTTTTTTTTTTT

     50     .    :    .    :    .    :    .    :    .    :
        TTTTTTTTTTTTTTTTTTTTTTAGTATTTCAGCAGGATCTGCTGGCAGGG
        ||||||||||||||||||||||  | |||        | |
        TTTTTTTTTTTTTTTTTTTTTT--T-TTT--------TTT----------

    100     .    :    .    :    .    :    .    :    .    :
        TTTTTTTGTTTTATTTGTTTGCTTATTTTTAAATTAACTGTTTTGAGCTT
        ||||||| |||| |||
        TTTTTTTTTTTTTTTT----------------------------------

    150
        TGA

        ---

21:3I2=3I1=3I2=1X40=4D1=7D7=1X4=1X3=1X3=
Before:Match:After:Start 45:39:66:15





# Scoring Matrix for edit distance:
"Ideally, the match/mismatch penalties used in genome alignment would match the
evolutionary distances of the sequences being aligned; human DNA to itself is
expected to be more than 99.9% identical."

"The match/mismatch ratios used in DNA similarity searches also have target
evolutionary distances. The stringent match/mismatch ratios used by MEGABLAST
are most effective at matching sequences that are essentially 100% identical,
e.g. mRNA sequences to genomic exons."
source:http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3848038/

RapMap should use a match/mismatch set to +1/-3 (roughly 99% similar genomes.)

Additional validation of edit distance:
"In general, the most widely used error models are the Hamming distance, which
accommodates only for mismatches between the read and a chosen genomic
location, and the edit distance, which accounts for mismatches and indels."
source:http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3627565/

# Alignment libraries
## [SSW](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/)
http://dx.plos.org/10.1371/journal.pone.0082138

Would need to configure to do semi-global alignment.

## [SeqAn](https://github.com/seqan/seqan/)
https://www.seqan.de/

Looks very promising but performance was prohibitive.


                              txpHitStart
                                  |
                                  v
TR: ACTTTTGTAAAGATTAAGCTCATTTAGTGTTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAGTATTTCAGCAGGATCTGCTGGCAGGGTTTTTTTGTTTTATTTGTTTGCTTATTTTTAAATTAACTGTTTTGAGCTTTGA

RC: GCCTTCCGTGCCTTGTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
FW: AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACACAAGGCACGGAAGGC
                   ^
                   |
                queryPos



Looks like the qa.pos is not valid???????????
Yes, in two places in the code we forgot to subtract the k-mer position from
the hit position. Instead of qa.pos pointing to the left most index of the read
into the transcript it pointed to the start of the match into the transcript.
This took some time to find because it was only when the read has only *one*
hit in the SA.



# Results

##### Single-end reads with no output (compiled with -O4, average over 5 runs)

| Threads  | Reads    | Quasimap  | Quasimap+Align | Slowdown |
|:--------:|:--------:| :--------:| :------------: | :------: |
|  1       |  1M      |   8.09s   |      31.47s    |   3.89x  |
|  8       |  1M      |   1.78s   |       6.60s    |   3.72x  |
|  1       |  10M     |   90.03s  |     369.99s    |   4.11x  |
|  8       |  10M     |   19.44s  |     84.51s     |   4.35x  |

##### Paired-end reads with no output (compiled with -O4, average over 5 runs)

| Threads  | Reads    | Quasimap  | Quasimap+Align | Slowdown |
|:--------:|:--------:| :--------:| :------------: | :------: |
|  1       |  1M      |   16.17s  |      58.11s    |   3.59x  |
|  8       |  1M      |    3.61s  |      12.73s    |   3.52x  |
|  1       |  10M     |  179.51s  |     701.82s    |   3.91x  |
|  8       |  10M     |   41.71s  |     164.21s    |   3.93x  |

# Improvements in RapMap Alignment

# Improvements in RapMap

Only reverse complement the read once at the start of the read. This could reduce code duplication and improve both readability and performance.
