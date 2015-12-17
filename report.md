# Efficiently Aligning RapMap Read Mappings
###### Shane Harvey (shane.harvey@stonybrook.edu)
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

### Quasi-mapping procedure
To understand the alignment procedure used later, one needs to understand the
basic quasi-mapping algorithm. Here is a simplified explanation, curious or
confused readers should refer to the pre-print for a detailed explanation and
illustration.

First, a read is scanned from left to right until a k-mer ki is encountered
that appears in the transcriptome (via the hash table). The hash table returns
the complete suffix array interval I(ki) where the ki is a prefix. This represents 
a suffix array hit.

After a k-mer hit we skip forward in the read to avoid redundant hits. The skip
is calculated as follows. The maximum mappable prefix (MMP) of the read to this
interval is computed. Finally, the next informative position (NIP) is computed
from the longest common prefix of the MMP interval. The search continues at the
position i + NIP - k in the read.

For our purposes it is important to note that the MMP of a hit a region of the
read that **exactly** matches the reference transcripts.

## Alignment using quasing-mapping
A naive approach would take the quasi-mapping locations for a read and align it
to the transcriptome region around each location. But we can leverage information
calculated during the quasi-mapping procedure to reduce or even eliminate the 
amount of the read alignment necessary.

Specifically we save the MMP and the start position of the k-mer. With this we
can calculate the region of the read and the transcript that match *exactly.*
This length of this region, the MMP, will be at least the size of the k-mer, but
in practice is often much larger. In our sample data with 76 base pair reads and
an index built with k-mers of size 31 the average MMP is ~70.

Below is a diagram of a sample read mapping with the exact match region shown.
```
          alignBefore                MMP                      alignAfter
               |                      |                            |
            |--v--| |-----------------v-------------------| |------v---------|
Transcript: GGCAGCC:GTGAGGCGCGTGTTCGGGCTCTTGCCGTCCCCGCACCCG:CACCGCGGTTACTGGCTT
Read:         ACGTA:GTGAGGCGCGTGTTCGGGCTCTTGCCGTCCCCGCACCCG:TTTTTT
                   ^
                   |
              k-mer positon
              
And the resultant optimal alignment:
           alignBefore              MMP                 alignAfter
                |                    |                       |
            |---v---|----------------v--------------------|--v--|
Transcript: GGCA-GCCGTGAGGCGCGTGTTCGGGCTCTTGCCGTCCCCGCACCCG------
               | |  |||||||||||||||||||||||||||||||||||||||      
Read:       ---ACGTAGTGAGGCGCGTGTTCGGGCTCTTGCCGTCCCCGCACCCGTTTTTT

And finally the equivalent CIGAR string:
1=1I1=2X39=6I
```

## Aligning before and after the match
To align the segment of the read before and after the match we implement semi-global
affine gap alignment. This was chosen because available libraries required
significant modifications or were too heavy weight.

##### SSW https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/
This is a very fast and highly optimized aligner but we would need to patch it
to do semi-global alignment and not Smith-Waterman. It would also
need to distinguish between match and mismatches in the traceback to report an
accurate CIGAR string.

##### SeqAn https://github.com/seqan/seqan/
SeqAn was initally used to perform alignments and output the CIGAR strings. It is
a very nice library and has excellent documentation. Although performance was
prohibitive and thus the need to implement alignment within RapMap, it was
extremely useful as reference implementation.

## Alignment scoring parameters
Research shows that the choice of parameters for the affine gap semi-global
alignment is very important. In the context where we are mapping RNA-seq data to
transcriptomes we should expect the transcripts and reads to be from very similar
genomes. As in MEGABLAST, we use match/mismatch scores of +1/-3 and like the SSW
library we use gap start/gap extend scores of -3/-1.

## Alignment results
RapMap now has support to quasi-align reads to its set of mapping locations.
This can be controlled with a command line flag, *-a*, which defaults to false.
If the flag is not present alignment is skipped and the before/after match segements
of the read are reported as **M** in the CIGAR string.

## Speed on synthetic data
To measure the speed of alignment we compare mapping speed between quasi-mapping with
and without alignment. As in the RapMap pre-print on bioRxiv,
we measure speed using a synthetic data set generated from the human
transcriptome. However due to hardware limitations, i.e. old laptops, we tested
using datasets of 1 million and 10 million 76 base pair, paired-end reads
instead of the ~48 million in the paper. The reference transcriptome consists
of 86,090 transcripts corresponding to protein-coding genes (same as the paper). 

Here single-end reads refer to looking at **only** the left-mates of the
paired-end reads. This was done to verify that single-end and paired-end
reads can be aligned properly.

On average we can expect a speed slowdown of 3.5-4.5x when running RapMap with
alignment verses quasi-mapping alone. Even still RapMap is faster than both
STAR and Bowtie 2. Also note that there are still significant optimizations
that can be done which will be discussed later.

All benchmarks were performed on a Intel(R) Core(TM) i7-2630QM CPU @ 2.00GHz with 8GB of RAM.
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

## Alignment improvements
There are two distinct ways that quasi-alignment can be substantially improved.
The first optimization would target the alignment DP itself. My implementation
of affine gap alignment leaves significant room for optimization. The
second route is to avoid repeating the exact same alignment. Both improvements
are described below.

### Avoid repeated alignments
We could remember the next informative position (NIP) of the match and use this information to avoid
repeating the same alignment after the match multiple times for a single read. Also this is complicated because the QuasiAlignment hits
could be for different parts of the SA ie they are hits from a different k-mer.
This means that care would need to be taken so that we don't compare the NIP from
hits originating from different k-mers in the read. This still ignores the part of
the transcript that matches with the read before the k-mer.

QuasiAlignment hits would need to be grouped by the originating SA interval hit
(could use the saved query position or *queryPos.*) Then we can see if the after
match alignment can be reused for each alignment in this group buy comparing the
saved NIP. This extra management should not be more complex than doing an extra
semi-global alignment.

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

Also it might
be useful to use a library such as SeqAn to manage SAM output. This would reduce
the scope of code RapMap needs to deal with by abstracting away output. Even
better, SeqAn has SAM and *BAM* support so we could output the user's preferred
format configured with a command line flag.

## RapMap bugs
During the addition of alignment serveral bugs were found and/or fixed. In one RapMap
consumed memory during *single-end* read mapping by not flushing the quasi-map hit buffer.
In another the reported hit position of the read was not correctly adjusted by the
length before the match.

## Conclusion
We present an extensiton to the tool RapMap that output the alignments of the reads to
the mapping locations. This degrades run time by 3.5-4.5x even still performing very
fast. Optimiaztions could reduce the overhead of alignment by avoiding repeats and
imporving the alignment algorithm.

## References
Avi Srivastava, Hirak Sarkar, Rob Patro (2015). RapMap: A Rapid, Sensitive and Accurate Tool for Mapping RNA-seq Reads to Transcriptomes. bioRxiv. http://dx.doi.org/10.1101/029652

Zhao M, Lee W-P, Garrison EP, Marth GT (2013) SSW Library: An SIMD Smith-Waterman C/C++ Library for Use in Genomic Applications. PLoS ONE 8(12): e82138. doi:10.1371/journal.pone.0082138

Andreas Döring, David Weese, Tobias Rausch and Knut Reinert. SeqAn an efficient, generic C++ library for sequence analysis. BMC Bioinformatics, 9:11, 2008.

Pearson, W. R. (2013). Selecting the Right Similarity-Scoring Matrix. Current Protocols in Bioinformatics / Editoral Board, Andreas D. Baxevanis ... [et Al.], 43, 3.5.1–3.5.9. http://doi.org/10.1002/0471250953.bi0305s43

Siragusa, E., Weese, D., & Reinert, K. (2013). Fast and accurate read mapping with approximate seeds and multiple backtracking. Nucleic Acids Research, 41(7), e78. http://doi.org/10.1093/nar/gkt005
