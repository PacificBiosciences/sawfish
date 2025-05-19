# Sawfish methods

## Candidate discovery
The first ‘discover’ step of sawfish is run once on each input sample. In this step each sample is scanned for SV and
CNV evidence. The SV evidence is clustered and used to assemble a set of candidate SV alleles. Details of this step are
provided below.

### Scanning sample alignments
The method scans over all mapped reads in the genome to identify clusters of breakpoint-associated read signatures and
to create regional depth bins. Reads are filtered from this scan if they are flagged as unmapped, secondary, QC-failed
or duplicate. Reads with gap-compressed identity less than 0.97 are filtered out as well. Reads with mapping quality
less than 10 are disqualified from breakpoint evidence scanning but still used to find regional sequencing depth.

Breakpoint evidence is gathered from indel and split-read annotations in the read alignments. In each case a simple
breakpoint candidate is created with breakends matching the location and orientation implied by the corresponding read
alignment feature. Split-read evidence is parsed from primary alignments only, with all supplementary alignment detail
taken from the auxiliary `SA` tag. Additional unpaired breakend evidence is gathered from soft-clipped read edges. Such
unpaired breakends are used to assemble large insertions when candidate breakend pairs are found in the expected
orientation (details below). Soft-clipped read ends only contribute to breakend candidates when at least 500 bases of
the read are soft-clipped on one end, and no clipping is found on the other end of the read.

Depth is enumerated from the alignments by recording the average depth within non-overlapping 1kb bins covering each
chromosome. Gaps created by splitting reads into primary and supplementary alignments are accounted for in the depth
calculation, as well as alignment deletions of 1kb or larger.

### Clustering breakpoint evidence
Breakpoint evidence from individual reads is clustered into candidate breakpoint clusters as follows: Each breakpoint
evidence observation is treated as a cluster with a supporting read count of 1. Breakpoint clusters with a matching
breakend orientation are tested for their total breakpoint distance, defined as the sum of the distance between each of
their breakends. If the breakpoint distance between the clusters is 500 bases or less, the breakpoint clusters are
merged, such that each merged breakend extends from the minimum of the two breakend start positions, and the maximum of
the two breakend end positions. The merged supporting read evidence of a cluster is the sum of the two input clusters.
After merging is completed, clusters with only a single supporting read are discarded as noise. All others comprise the
candidate breakpoint cluster set.

### Candidate breakpoint cluster refinement
In the breakpoint cluster refinement process, candidate breakpoints are assembled into SV haplotype contigs, which are
then aligned back to the genome to generate candidate SVs used in downstream merging and genotyping steps.

Refinement begins by defining the regions of the genome used to extract reads for SV contig assembly for each breakpoint
candidate. For breakpoint candidates with distant breakends, these will simply be the two breakend regions. For
breakpoint candidates that form an indel-like breakend orientation pattern with a breakend distance of 600 bases or
less, the assembly region is merged to span the full region between the two breakends. These become single-region
candidate breakpoints. Such breakpoints go through an additional clustering step to consolidate all single-region
candidate breakpoint regions within 300 bases of each other into a single candidate assembly region, except that the
consolidation process is limited to prevent the creation of consolidated regions larger than 8000 bases.

### Large insertion candidate generation
Large insertions can only be discovered by the standard breakpoint clustering process if the read mapper represents
large insertions in the reported read alignments. Additional candidate large insertions can be identified from local
soft-clipped read alignment patterns as follows. When a pair of left-anchored and right-anchored soft-clipped breakend
candidates are found such that the left and right breakends are within 500 bases, these are converted into large
insertion candidate assembly regions if they aren’t already overlapping a candidate assembly region from the standard
cluster refinement process.

### Consensus contig generation
The first step of contig generation is obtaining the subsequences of the reads around each candidate breakpoint
position. To do so, reads mapped near the candidate breakpoint are enumerated with the same filtration used for
breakpoint evidence discovery. For each read, a trimmed subread is extracted comprising 300 bases of read sequence in
each direction from the last base mapped on each side of the candidate breakpoint location. These trimmed reads are
additionally filtered out if they lack breakpoint evidence, requiring either an indel at least 25 bases long, or a
soft-clipped segment within the trimmed region. Note for the special case of large insertion candidates discovered from
paired breakend clusters, no read trimming is applied on the soft-clipped side of the candidate breakend, because the
length of the large insertion isn’t known for this case. If there are more than 100 trimmed reads for a given breakpoint
candidate, the reads are deterministically subsampled to 100.

The trimmed reads are next clustered and polished into a set of contigs using a simple iterative procedure to assign
each read to a partial order alignment (POA) graph using the spoa library (https://github.com/rvaser/spoa). Each of
these graphs are generated from previously evaluated trimmed reads, and each is interpreted as being sampled from the
same underlying haplotype combined with sequencing noise. Each trimmed read is successively aligned to each POA graph
using a linear gap alignment with weights of 1, -3, -1 for match, mismatch, and gap, using wavefront alignment
(https://github.com/smarco/WFA2-lib). If at least one read alignment has an aligned read length of 100 or higher, and an
alignment score normalized by the aligned read length of 0.96 or higher, then the trimmed read is assigned to the POA
graph with the highest normalized alignment score, and the POA graph is updated to include the new trimmed read. If the
read does not have a sufficiently high-quality alignment to any existing POA graph, then the read is used to start a new
POA graph so long as this would not create more than 8 graphs, otherwise the read is filtered out of the consensus
contig generation procedure. After all trimmed reads have been processed, the POA graphs are filtered to remove cases
with less than 2 supporting reads, and the top P remaining POA graph clusters by supporting read count are used to
generate a consensus contig for downstream steps, where P is the local ploidy count.

### Candidate SV generation
Top contigs generated for each breakpoint cluster are processed for candidate SV information by locally aligning these
back to the expected reference locations. For indel-like contigs assembled from reads over a single region of the
reference, the contigs can simply be aligned back to a similar segment of the reference sequence and all indels over the
minimum size (35 by default) extracted from the alignment as SV candidates.

For all other types of contigs, a synthetic derived chromosome reference is created by appending the two reference
region segments corresponding to each candidate breakend location, possibly with one of the reference segments
reverse-complemented according to the candidate breakend orientation pattern. In this way the contig can be aligned to
the synthetic reference with an apparent large deletion extending between the two reference regions, and this alignment
can be processed back into a pair of breakends of any orientation, which generalizes to handle all types of SV
breakpoints.

Whether using the small indel or generalized breakpoint reference alignment procedure, each breakpoint alignment can be
standardized by left-shifting its position and finding the full breakend homology range and insertion sequence.

### Depth segmentation
In addition to SV candidate generation, sawfish also performs depth segmentation during the discover step. This step is
intended to estimate a 'draft' of the sample segmentation, which will be used to estimate the haploid coverage of the
sample. This 'draft' segmentation will be further refined in the joint-call step prior to integration with
breakpoint-based SV calls and final variant output.

The first step of the depth segmentation procedure is to estimate GC-bias in the sample so that this confounder can be
accounted for in the segmentation procedure. GC-bias estimation is followed by an iterative segmentation procedure where
the method alternates between generating a copy-number segmentation for the sample and then using that segmentation to
estimate haploid coverage. This iteration is continued until the haploid coverage estimate stabilizes. While this
procedure has little impact on the haploid coverage estimate of typical samples, it provides an important correction for
samples with very large copy variant regions.

#### GC-bias correction
The GC content of each depth bin is taken from the reference sequence, and used to implement a simple correction factor
in the segmentation process as described below.

For each depth bin, the GC and AT counts are taken from the reference sequence within a 20kb window centered on the
depth bin. The GC counts window size is independent of the depth bin size. It is adjustable with hidden option
`gc-genome-window-size`. GC counts are also accumulated for the entire genome.

Next the method iterates over all depth bins to build a mapping between GC fraction and average depth. For this
operation the GC fraction is discretized into equal size frequency bins, using 40 bins by default (adjustable with the
option `gc-level-count`). Each depth bin is assigned a GC level, and the depth for all bins at a given GC level is
accumulated and averaged. Any GC-level with at least 2 million bases of depth bin representation is qualified to be used
in GC bias correction process. For example, for a depth bin size of 1000, a given GC level would need at least 2000
depth bins observed at that GC level to be qualified for the GC bias correction process. The qualifying GC level with
the highest average depth is found, and the depth at this GC-level is interpreted as being the 'unbiased' coverage
level, and a correction factor is computed for each GC level by taking the average depth for the GC level and dividing
it by the 'unbiased' coverage level from the max depth bin.

For non-qualifying GC-level bins, without sufficient depth bin support, the  correction factor is found as follows. Find
the range of continuously qualifying GC-level bins containing the max-depth GC-level. All bins higher than that range
use the correction factor of the highest qualifying bin in the range, all bins lower than that range use the correction
factor of the lowest qualifying bin in the range. The procedure for computing correction factors is implemented in
`get_gc_correction`.

GC bias correction is applied downstream emission probabilities in the depth segmentation model. the depth values
themselves are never adjusted.

#### Segmentation method
Sample haploid depth is estimated from a subset of chromosomes which are assumed to be diploid and represent
high-quality assemblies. By default this subset is found from chromosomes names which match a regular expression
intended to match human autosomes. The haploid depth is computed from the zero-excluded mean depth of this chromosome
set. When GC content correction is applied, the haploid depth is adjusted to reflect the haploid coverage level of the
least-biased GC window. When available, the depth estimation process will use copy number state estimates from a
previous round of segmentation to improve accuracy.

Copy number segmentation is run using a hidden markov model (HMM). The allowed copy number states include copy numbers 0
to 4 and a 'High' state intended to model copy numbers 5 and higher, there is also a special 'Unknown' copy number state
used to implement the method's excluded region behavior.

The HMM is setup by default with a single transition probability of 0.02 applied for any change between states. In
subsequent SV/CNV merging steps, the transition probabilities for certain bins and copy number state change directions
(increasing or decreasing) will be increased to reward copy number transitions at putative SV breakend locations.

The emission probability of the standard copy-number states is a Poisson sampling of the observed depth given an
expected depth equal to the haploid coverage estimate times the state's associated copy number, divided by the estimated
GC-bias.

For depth bins in excluded regions, the emission probability of the 'Unknown' state is 1 and and of any other state is
0.5 (prior to normalization), for depth bins in non-excluded regions the emission probability of the 'Unknown' state is
0. This scheme allows for a continuous copy number segment to span a small excluded region, while producing 'Unknown'
state segments across larger excluded regions.

A modification to the state emission probabilities is needed to account for the depth bin size being smaller than the
average read size, resulting in a high degree of dependency between the depth observations at adjacent bins, violating
standard HMM assumptions. A modeled correction with conditional emission probabilities has been explored for the purpose
of correcting for this problem via conditional independence, but in practice we find a simple empirical adjustment
is effective: all emission probabilities are raised to the power of the dependency correction factor, which is
0.02.

Given the above transition and emission structure, segmentation is performed by a Viterbi parse of the depth bins. The
segmentation process does not consider any small-variant minor allele frequency input.

During the discover step, the combined process of haploid depth estimation and segmentation is iterated until the depth
estimate converges.

## Joint calling
The second ‘joint-call’ step of the sawfish pipeline enables candidate SVs and CNVs to be analyzed across multiple
samples. The primary steps are:

1. Consolidate SV candidates which are duplicated across samples
2. Evaluate read support for the deduplicated candidate SV haplotypes, to genotype SVs across all samples
3. Merge the joint-genotyped SV set with copy number segmentation results
    - Reward copy number transitions at selected SV breakpoints
    - Re-segment depth in all samples with breakpoint-adjusted transition probabilities
    - Segmentation 'sync' step: reward small adjustments in copy number boundary locations to match across samples
    - Merge matching SV/CNV events across all samples
    - Filter large unbalanced SVs without CNV support into breakpoints.
4. Report all SV and CNV results for all samples to a single integrated VCF output

### Duplicate haplotype consolidation
As a first step to duplicate haplotype merging, overlapping haplotypes from all samples are pooled into groups from
which duplicates are found and consolidated. This pooling procedure is run separately within the set of indel-like SV
candidates consolidated to a single reference region (as described in the candidate discovery section above), and all
other candidates associated with multiple reference regions. For indel-like candidates, the candidate pool is found from
all intersecting candidate regions. For all other SV candidates associated with two reference regions, pools are
composed of candidates where both reference regions are within 100 bases of at least one other candidate in the pool.

Within each candidate haplotype pool, candidates are clustered into duplicate groups based on pairwise testing first for
matching breakpoints, then for very high haplotype sequence similarity if the breakpoints aren’t an exact match. An
exception is made to exclude duplicate merging of haplotypes candidates from the same sample. Haplotype sequence
similarity is determined by using a linear gap aligner with 1, -3, -2 for match, mismatch, and gap scores. If the
resulting alignment score normalized by the aligned haplotype length is at least 0.97, then the haplotypes are treated
as duplicates. Within each duplicate haplotype group, one member is chosen as representative based on having the highest
supporting read count from contig assembly, or longest contig size when supporting read count is tied.

### Determining supporting read counts for each allele
SV allele read support is determined in the context of the overlapping haplotype pools created for the purpose of
duplicate haplotype identification. Each breakend of each SV allele is evaluated by aligning segments of locally mapped
reads to the corresponding segment of the SV haplotype assembly, in addition to aligning these reads to the reference
sequence and other SV alleles in the overlapping haplotype pool. The segments of the read, reference and SV haplotypes
selected for this purpose on each candidate SV breakend are extended up to 500 bases from the non-anchored edge of the
breakend’s homology range. Alignment is scored using a linear gap aligner with 1, -3, -2 for match, mismatch, and gap
scores, with alignment scores normalized by the aligned read segment length.

The alignment scores for each read are next converted into support counts, where each read could be identified as
uniquely consistent with either the haplotype of the SV allele in question, the reference haplotype or a candidate
overlapping SV haplotype. Candidate overlapping haplotypes are only identified for local indel-like SV candidates. The
overlapping case can be from a second SV haplotype assembled from reads in the given sample. If a second haplotype
wasn’t assembled for the sample, then overlapping SV haplotypes from other samples (among those remaining after merging
duplicate haplotypes) are compared to find the overlapping SV haplotype with the highest support from all sample reads.
If such an overlapping case is found it is added as a second ‘guest’ candidate haplotype.

Each read’s support for the SV haplotype is evaluated separately for each breakend of a given breakpoint. A read is
counted as supporting an SV allele if its alignment score to the SV haplotype is better than the reference haplotype for
any breakend. A read can only support the reference haplotype if the reference haplotype alignment has a higher score on
all evaluated breakends. This arrangement helps to counteract various forms of reference bias in the alignment scores.
Reads supporting an overlapping SV haplotype are internally recorded for the overlapping haplotype but converted into
reference allele support for both quality score calculations and in the final VCF allele count output, per the reporting
convention used by most available SV calling tools.

### Merging of SV candidates with copy number segmentation results

After finding the supporting read count for each SV breakpoint allele, this breakpoint-based SV information is
integrated with the depth-track segmentation process over several steps as described below.

#### Rewarding copy number transitions at selected SV breakpoints

An increased copy number state transition probability is set for depth bins corresponding to certain SV breakends, prior
to re-segmenting the depth in all samples. The breakends where these transition bonuses a provided are limited to those
meeting the criteria described below to prevent excessive fragmentation of the copy number segments in more difficult
regions of the genome.

Breakpoints are used to increase copy number transition probabilities if they would not already be marked as filtered
in the VCF output, and one of the following conditions are met:

1. The breakpoint's orientation is consistent with that of a deletion or duplication of at least 50kb, and its span
intersects a constant copy number segment which extends at least 5kb past the ends of the breakpoint span. Providing a
bonus in this circumstance can help detect relatively smaller CNVs, where the depth change may be too subtle to segment
using the default transition scores, but using the breakpoint in question as a prior allows for increased CNV
sensitivity.

2. The breakpoint's orientation is consistent with that of a deletion or duplication of at least 50kb, its span overlaps
a CNV segment with at least an 80% reciprocal overlap, the copy number change is consistent with the breakpoint genotype
and orientation, and only one such breakpoint meets these conditions for the CNV segment in question. For SV/CNV pairs
which are close to meeting the merge threshold for reciprocal overlap the copy number transition bonus provided by this
rule can help delineate whether the size difference is primarily explained by depth sampling noise of the same
underlying variant or by two different variants in the sample.

The bonus transition probability is the square root of the standard transition probability, restricted to only copy number
gains or losses across the SV breakend bin as would be expected from the breakend direction. After recording all the bonus
transition probability bins in every sample, then depth for all samples goes through a final re-segmentation step.

#### Multi-sample copy-number boundary refinement

To help identify when a highly similar CNV segment is shared across multiple samples, the multi-sample segmentation
result is subject to a simple multi-sample breakpoint refinement step, where each copy number boundary is allowed to
shift by up to 5kb to match a copy number boundary of the same orientation in another sample. The probabilistic penalty
of this shift is computed from the emission score differences over the shifted copy-number boundary bins, and if the
shifting cost is less than a small synchronization bonus then the shift is retained.

The idea of applying a small synchronization bonus is to bring copy number boundaries together when the same underlying
event across multiple samples is represented with different boundaries due to sampling noise, rather than cases
representing different variants.

#### Merging SV and CNVs across multiple samples

After re-segmenting depth in all samples and allowing copy number boundaries to have a small synchronization shift, the
final merging process determines which SV and CNV events will be merged together in the final VCF output -- this process
will also determine which large unbalanced SVs will be converted to breakends due to the lack of a matching CNV. SVs can
be merged with up to one CNV per sample and must be the only matching SV for each given CNV. The criteria for the SV and
CNVs from each sample to merge are as follows:

- SVs are eligible for merging with CNVs if the SV's breakpoint pattern is consistent with deletion or duplication, if
the breakpoint span is at least 10kb, and if the SV has not met the criteria for any output filtration flags.
- CNVs are eligible for merging with SVs if the CNV has not met the criteria for any output filtration flags.
- SV and CNV spans must have a reciprocal overlap of at least 90% to be merged.
- The edges of the SV and CNV spans can be separated by no more than 8kb to be merged. This maximum separation criteria
is required at both the left and right edges of the span.
- The SV and per-sample matching CNVs must have a compatible type in at least 50% of samples where the SV or CNV has a
variant genotype or copy number. Compatible type means that SV deletions are paired with a copy loss and SV duplications
are paired to a copy gain. Homozygous SV deletions must be more specifically paired to a zero copy number segment.
- Where multiple SVs are matched to the same CNV, the SVs are ranked based on the most samples in which the SV has a
candidate CNV match, and secondarily the least total SV/CNV span difference over all candidate CNV matches. Only the top
ranked SV is retained as eligible for merging with the CNV.

Any remaining SV-CNV pairings after applying all of these criteria are merged and reported in a unified VCF record
format.

#### Large SV filtration to breakpoints

Any SV deletions and duplications spanning 50kb or larger that have not been paired to a CNV in any sample will be
reported in out as a pair of VCF breakend (BND) records.


#### Generate copy number quality scores

For each copy number segment delineated in one or more samples, we find the probability that the segment copy number
is equal to the segmented copy number.

### Quality models

Sawfish includes a standard quality model applied to all structural variants based on read support for the SV
breakpoint. These qualities are the basis for the VCF `GQ` (per-sample) and `QUAL` (per-variant) values for all
breakpoint-based SVs.

An additional quality model is described below which expresses the strength of depth-support evidence for the copy
number called in each CNV, these qualities are the basis for the VCF `CNQ` (per-sample) confidence value for all
variants with depth-support(both merged and unmerged CNVs). For unmerged CNVs the depth support model is also used to
compute the `QUAL` value for the CNV overall.

For merged SV/CNV records, the `QUAL` value reflects the higher of the two (breakpoint-based or depth-based) quality
scores.

#### Read-support quality model
The sawfish SV quality scores are generated from a simplified model which generates qualities directly from the read
support counts for each allele. Ideally the scores of read alignments to each allele haplotype would be used instead of
counts, but thus far the count approximation hasn’t been attributed to any substantial fraction of SV genotyping errors.

For each SV allele and each sample, we solve for the posterior probability of diploid genotypes, $G$, given the observed
read counts at the SV loci, $D$ as follows

$$
P(G \vert D) \propto P(D \vert G) P(G)
$$

The diploid genotypes are $G = \{\text{ref}, \text{het}, \text{hom}\}$ representing 0, 1 or 2 copies of the SV allele.
The genotype prior is

$$
P(G)=
\left\lbrace
  \begin{array}{rl}
\theta      & \mbox{ if het} \\
\theta / 2  & \mbox{ if hom} \\
1 - \theta 3 / 2  & \mbox{ if ref}
\end{array}
\right.
$$

where $\theta = 5 \times 10^{-4}$.

The genotype likelihood is found as follows

$$
P(D \vert G) = \prod_{d \in D} P(d \vert G)
$$

treating each read observation $d \in D$ as independent. The read likelihood is

$$
P(d \vert G) = \sum_{a \in A} P(d \vert a) P(a|G)
$$

where $A = \lbrace\text{ref}, \text{alt}\rbrace$ are the SV alleles in the model representing the reference and SV
haplotypes. As previously discussed, support for any overlapping SV haplotypes is counted towards the reference allele
for the purpose of the genotype quality model, which is intentionally simplified to represent only one SV haplotype at a
time.

Considering the terms in the read likelihood, the allele likelihood $P(d \vert a)$ is set from the read allele support
counts using a single erroneous read support count probability $e = 1 \times 10^{-5}$ for all cases. The allele
probabilities $P(a \vert G)$ are the simple allele fractions (0, 0.5, 1) associated with each genotype.

#### Depth-support quality model
The quality model used for CNVs is based on the probability of different copy number values, given a particular segment
and haploid depth estimate. Conditioning CNV qualities on fixed segmentation and haploid depth values are simplifying
assumptions that should be considered when interpreting these qualities.

Similar to the read-based quality model, for the copy number model we seek to approximate the posterior over copy number
states $C$ given the depth information from the sequencing reads.

$$
P(C \vert D) \propto P(D \vert C) P(C)
$$

For simplicity a naive prior is used for the copy number states.

The copy number likelihood is

$$
P(D \vert C) = \prod_{d \in D} P(d \vert C)
$$

where $d \in D$ represents the depth bins across the full CNV segment. The likelihood for each depth bin is the same as
the dependency corrected HMM emission probability discussed above in the context of depth segmentation. This is

$$
P(d \vert C) = \text{PoissonPMF}(c_{\text{obs}};\lambda = c_{\text{expect}})^{x_{\text{dep}}}
$$

Here $c_{\text{obs}}$ is observed coverage in depth bin $d$ and the expected coverage is $c_{\text{expect}} =
c_{\text{hap}} n g$, given $c_{\text{hap}}$ is the haploid depth estimate, $n$ is the copy number associated with each
copy number state in $C$, and $g$ is the depth reduction factor estimated due to GC-bias for depth bin $d$. The
additional term $x_{\text{dep}}$ is the dependency correction factor. This was introduced to approximately correct for
the dependency between depth bins due to the read length being substantially longer than the depth bin.
