# Efficiently Aligning RapMap Read Mappings
RapMap quickly reports hit locations of sequencing reads to a transcriptome. In
this report we extend RapMap to efficiently compute the optimal alignment of the
read to the transcriptome at these hit locations. As a result we only run
3.5-4.5x slower to generate the alignment CIGAR strings. The speed comes from
leveraging the exact match region found in quasi-map lookup and significant
optimizations are still possible.

## RapMap
RapMap is an implementation of quasi-mapping, in which sequencing reads are
mapped to a target transcriptome. Quasi-mapping is a novel algorithm which uses
a generalized suffix array (SA) and a hash table to efficiently and accurately
determine the likely origin locations of the sequencing read. Using quasi-mapping
RapMap is capable of mapping reads considerably faster than existing tools.

### Quasi-mapping procedure [TODO Simplify]
To understand the alignment procedure used later one needs to understand the
basic quasi-mapping algorithm. Here is a simplified explanation, curious or
confused readers should refer to the pre-print for a detailed explanation and
illustration.

First, a read is scanned from left to right until a k-mer ki is encountered
that appears in the transcriptome (via the hash table.) The hash table returns
the complete suffix array interval I(ki) where the ki is a prefix. This represents 
a suffix array hit.

After a k-mer hit we skip forward in the read to avoid redundant hits. The skip
is calculated as follows. The maximum mappable prefix (MMP) of the read to this
interval is computed. Finally, the next informative position (NIP) is computed
from the longest common prefix of the MMP interval. The search continues at the
position i + NIP - k in the read.

For our purposes it is important to note that the MMP of a hit a region of the
read that **exactly** matches the reference transcripts.

## Alignment Intro [TODO]
A naive approach would take the quasi-mapping locations for a read and align it
to the transcriptome region around each location. But we can leverage information
calculated during the quasi-mapping procedure to reduce or even eliminate the 
amount of the read alignment necessary.

Specifically we save the MMP and the start position of the k-mer. With this we
can calculate the region of the read and the transcript that match *exactly.*
This length of this region, the MMP, will be at least the size of the k-mer, but
in practice is often much larger. In our sample data with 76 base pair reads and
an index built with k-mers of size 31 the average MMP in ~70.

Below is a diagram of a sample read mapping with the exact match region shown.


        txpHitStart
          |
          |txpStart                   txpHitEnd             txpEnd
          |   |txpMatchStart  txpMatchEnd |                    |
          v   v    v               v      v                    v
    ----------|----================----------------------------|---------|
          |---|----================|------|
                   ^               ^
                   |            matchEnd
               matchStart


## Aligning before and after the match
To align the segment of the read before and after the match we implement semi-global
affine gap alignment. This was chosen because available libraries required
significant modifications (SSW) or were too heavy weight (SeqAn.)

## Scoring matrix for edit distance:
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

## Alignment results
RapMap now has support to quasi-align reads to its set of mapping locations.
This can be controlled with a commandline flag, *-a*, which defaults to false.
If the flag is not present Gotoh alignment is skipped and the would be aligned
suffix and/or prefix of the read are reported as **M** in the CIGAR string.

## Speed on synthetic data
To measure the of alignment we compare mapping speed between quasi-mapping with
and without alignment. As in the [RapMap pre-print on bioRxiv](http://biorxiv.org/content/early/2015/10/28/029652),
we measure speed using a synthetic data set generated from the human
transcriptome. However due to hardware limitations, i.e. old laptops, we tested
using datasets of 1 million and 10 million 76 base pair, paired-end reads
instead of the ~48 million in the paper. The reference transcriptome consists
of 86,090 transcripts corresponding to protein-coding genes (same as the paper.) 

Here single-end reads refer to looking at **only** the left-mates of the
paired-end reads. This was done to verify that single-end and paired-end
reads can be aligned properly.

On average we can expect a speed slowdown of 3.5-4.5x when running RapMap with
alignment verses quasi-mapping alone. Even still RapMap is faster than both
STAR an Bowtie 2. Also note that there are still significant optimizations
that can be done which will be discussed later.

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

## RapMap bugs [TODO]
SAM flags field is wrong (0x900 for only first read mapping)? Perhaps it would
be nice to use a library such as SeqAn to manage SAM output. This would reduce
the scope of code RapMap needs to deal with by abstracting away output. Even
better SeqAn has SAM and *BAM* support so we could output the user's preferred
format configured with a command line flag.

```
                              txpHitStart
                                  |
                                  v
TR: ACTTTTGTAAAGATTAAGCTCATTTAGTGTTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAGTATTTCAGCAGGATCTGCTGGCAGGGTTTTTTTGTTTTATTTGTTTGCTTATTTTTAAATTAACTGTTTTGAGCTTTGA

RC: GCCTTCCGTGCCTTGTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
FW: AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACACAAGGCACGGAAGGC
                   ^
                   |
                queryPos
```

Looks like the qa.pos is not valid???????????
Yes, in two places in the code we forgot to subtract the k-mer position from
the hit position. Instead of qa.pos pointing to the left most index of the read
into the transcript it pointed to the start of the match into the transcript.
This took some time to find because it was only when the read has only *one*
hit in the SA.

## Alignment improvements
There are two distinct ways that quasi-alignment can be substantially improved.
The first optimization would target the alignment DP itself. My implementation
of Gotoh affine gap alignment leaves significant room for optimization. The
second route is to avoid repeating the exact same alignment. Both improvements
are described below.

### Optimize Gotoh affine gap alignment [TODO]
Group matrices (m, x, y, tm, tx, ty) for better cache performance?

#### Alignment libraries
##### [SSW](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/)
http://dx.plos.org/10.1371/journal.pone.0082138

Would need to configure to do semi-global alignment and not Smith-Waterman. Also
need to distinguish between match and mismatches in the traceback to report an
accurate CIGAR string.

##### [SeqAn](https://github.com/seqan/seqan/)
https://www.seqan.de/

Looks very promising but performance was prohibitive. 

### Avoid repeated alignments [TODO]
Cache transcript/read to avoid repeated alignments! This could lead to a big
performance increase because often the hit regions are exactly the same bewteen
transcripts.

Or could remember the Next Informative Position of the match (but what abut the
before the match part?) Also this is complicated because the QuasiAlignment hits
could be for different parts of the SA ie they are hits from a different k-mer.
This means that care would need to be taken so that we don't compare the NIP from
hits originating from different k-mers in the read. This still ignores the part of
the transcript that matches with the read before the k-mer.

QuasiAlignment hits would need to be grouped by the originating SA interval hit
(could use the saved query position or *queryPos.*) Then we can see if the after
match alignment can be reused for each alignment in this group buy comparing the
saved NIP. This extra management should not be more complex than doing an extra
semi-global alignment.

Test cases to make sure aligner and CIGAR reporting works properly.


## RapMap improvements
There are a number of improvements to be made in RapMap both in terms of
readability and performance. First, we could only reverse complement the read
once at the start of processing a read. This would improve performance because
a read would not need to be reverse complemented in multiple places leading to
wasted memory allocation and CPU cycles.

This change could additionally reduce code duplication and improve readability, performance, and
maintainability. Readability because there would less code to read. Performance
because a read would not need to be reverse complemented in multiple places.
Maintainability because we would have less chance of modifying the code in copy
A but forgetting to update copy B.

```c++
void collectHits(string& read, vector<QuasiAlignment>& hits, /* params */) {
      /* A) code to process read k-mers */
      ...
      /* B) code to process reverse complement read k-mers */
      ...
      /* A) code to process read hits */
      ...
      /* B) code to process reverse complement read hits */
      ...
}
```
An refactored version would look something like,
```c++
void collectHits(string& read, vector<QuasiAlignment>& hits, /* params */) {
      /* code to process read k-mers */
      ...
      /* code to process read hits */
      ...
}

void processRead(string& read) {
      string revCompRead = reverseComplement(read);
      vector<QuasiAlignment> hits;
      collectHits(read, hits, /* params */);
      collectHits(revCompRead, hits, /* params */);
      
      /* output hits */
}
```

This pattern is present in multiple places in the RapMap code base. These
proposed changes would improve code health making it easier for fresh minds
to become familiar with the RapMap code base. Additionally it would be easier
to contribute new features and optimizations, because of de-duplication.

## Future Work [TODO]

## Conclusion [TODO]

## Some alignments are sub-optimal! [TOOD verify]
```
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
```

```
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
```
