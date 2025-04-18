# Sawfish methods

## Candidate discovery
The first ‘discover’ step of sawfish is run once on each input sample. In this step each sample is scanned for SV evidence, which are clustered and used to assemble a set of candidate SV alleles. Details of this step are provided below.

### Scanning sample alignments
The method scans over all mapped reads in the genome to identify clusters of breakpoint-associated read signatures and to create regional depth bins. Reads are filtered from this scan if they are flagged as unmapped, secondary, QC failed or duplicate. Reads with gap compressed identity less than 0.97 are filtered out as well. Reads with mapping quality less than 10 are disqualified from breakpoint evidence scanning but still used to find regional sequencing depth.

Breakpoint evidence is gathered from indel and split read annotations in the read alignments. In each case a simple breakpoint candidate is created with breakends matching the location and orientation implied by the corresponding read alignment feature. Note that split read evidence is parsed from primary alignments only. Additional unpaired breakend evidence is gathered from soft-clipped read edges. Such unpaired breakends are used to assemble large insertions when candidate breakend pairs are found in the expected orientation (details below). Soft-clipped read ends only contribute to breakend candidates when at least 500 bases of the read are soft-clipped on one end, and no clipping is found on the other end of the read.

The average depth of each 2kb bin across the genome is also found while scanning alignments for breakpoint evidence. Gaps created by splitting the read into primary and supplementary alignments are accounted for in the depth calculation, but not the alignment indels.

### Clustering breakpoint evidence
Breakpoint evidence from individual reads is clustered into candidate breakpoint clusters as follows: Each breakpoint evidence observation is treated as a cluster with a supporting read count of 1. Breakpoint clusters with a matching breakend orientation are tested for their total breakpoint distance, defined as the sum of the distance between each of their breakends. If the breakpoint distance between the clusters is 500 bases or less, the breakpoint clusters are merged, such that each merged breakend extends from the minimum of the two breakend start positions, and the maximum of the two breakend end positions. The merged supporting read evidence of a cluster is the sum of the two input clusters. After merging is completed, clusters with only a single supporting read are discarded as noise. All others comprise the candidate breakpoint cluster set.

### Candidate breakpoint cluster refinement
In the breakpoint cluster refinement process, candidate breakpoints are assembled into SV haplotype contigs, which are then aligned back to the genome to generate candidate SVs used in downstream merging and genotyping steps.

Refinement begins by defining the regions of the genome used to extract reads for SV contig assembly for each breakpoint candidate. For breakpoint candidates with distant breakends, these will simply be the two breakend regions. For breakpoint candidates that form an indel-like breakend orientation pattern with a breakend distance of 600 bases or less, the assembly region is merged to span the full region between the two breakends. These become single-region candidate breakpoints. Such breakpoints go through an additional clustering step to consolidate all single-region candidate breakpoint regions within 300 bases of each other into a single candidate assembly region, except that the consolidation process is limited to prevent the creation of consolidated regions larger than 8000 bases.

### Large insertion candidate generation
Large insertions can only be discovered by the standard breakpoint clustering process if the read mapper represents large insertions in the reported read alignments. Additional candidate large insertions can be identified from local soft-clipped read alignment patterns as follows. When a pair of left-anchored and right-anchored soft-clipped breakend candidates are found such that the left and right breakends are within 500 bases, these are converted into large insertion candidate assembly regions if they aren’t already overlapping a candidate assembly region from the standard cluster refinement process.

### Consensus contig generation
The first step of contig generation is obtaining the subsequences of the reads around each candidate breakpoint position. To do so, reads mapped near the candidate breakpoint are enumerated with the same filtration used for breakpoint evidence discovery. For each read, a trimmed subread is extracted comprising 300 bases of read sequence in each direction from the last base mapped on each side of the candidate breakpoint location. These trimmed reads are additionally filtered out if they lack breakpoint evidence, requiring either an indel at least 25 bases long, or a soft-clipped segment within the trimmed region. Note for the special case of large insertion candidates discovered from paired breakend clusters, no read trimming is applied on the soft-clipped side of the candidate breakend, because the length of the large insertion isn’t known for this case. If there are more than 100 trimmed reads for a given breakpoint candidate, the reads are deterministically subsampled to 100.

The trimmed reads are next clustered and polished into a set of contigs using a simple iterative procedure to assign each read to a partial order alignment (POA) graph using the spoa library (https://github.com/rvaser/spoa). Each of these graphs are generated from previously evaluated trimmed reads, and each is interpreted as being sampled from the same underlying haplotype combined with sequencing noise. Each trimmed read is successively aligned to each POA graph using a linear gap alignment with weights of 1, -3, -1 for match, mismatch, and gap, using wavefront alignment (https://github.com/smarco/WFA2-lib). If at least one read alignment has an aligned read length of 100 or higher, and an alignment score normalized by the aligned read length of 0.96 or higher, then the trimmed read is assigned to the POA graph with the highest normalized alignment score, and the POA graph is updated to include the new trimmed read. If the read does not have a sufficiently high-quality alignment to any existing POA graph, then the read is used to start a new POA graph so long as this would not create more than 8 graphs, otherwise the read is filtered out of the consensus contig generation procedure. After all trimmed reads have been processed, the POA graphs are filtered to remove cases with less than 2 supporting reads, and the top P remaining POA graph clusters by supporting read count are used to generate a consensus contig for downstream steps, where P is the local ploidy count.

### Candidate SV generation
Top contigs generated for each breakpoint cluster are processed for candidate SV information by locally aligning these back to the expected reference locations. For indel-like contigs assembled from reads over a single region of the reference, the contigs can simply be aligned back to a similar segment of the reference sequence and all indels over the minimum size (35 by default) extracted from the alignment as SV candidates.

For all other types of contigs, a synthetic derived chromosome reference is created by appending the two reference region segments corresponding to each candidate breakend location, possibly with one of the reference segments reverse-complemented according to the candidate breakend orientation pattern. In this way the contig can be aligned to the synthetic reference with an apparent large deletion extending between the two reference regions, and this alignment can be processed back into a pair of breakends of any orientation, which generalizes to handle all types of SV breakpoints.

Whether using the small indel or generalized breakpoint reference alignment procedure, each breakpoint alignment can be standardized by left-shifting its position and finding the full breakend homology range and insertion sequence.

## Joint calling
The second ‘joint-call’ step of the sawfish pipeline enables candidate SVs to be analyzed across multiple samples. The primary steps are: (1) consolidation of SV candidates which are duplicated across samples; (2) evaluation of read support for the deduplicated candidate SV haplotypes to genotype the SV in each sample; (3) evaluation of depth support for large copy-changing events; and (4) reporting all SVs as a VCF jointly genotyped across all samples.

### Duplicate haplotype consolidation
As a first step to duplicate haplotype merging, overlapping haplotypes from all samples are pooled into groups from which duplicates are found and consolidated. This pooling procedure is run separately within the set of indel-like SV candidates consolidated to a single reference region (as described in the candidate discovery section above), and all other candidates associated with multiple reference regions. For indel-like candidates, the candidate pool is found from all intersecting candidate regions. For all other SV candidates associated with two reference regions, pools are composed of candidates where both reference regions are within 100 bases of at least one other candidate in the pool.

Within each candidate haplotype pool, candidates are clustered into duplicate groups based on pairwise testing first for matching breakpoints, then for very high haplotype sequence similarity if the breakpoints aren’t an exact match, with an exception made to exclude duplicate merging of haplotypes candidates from the same sample. Haplotype sequence similarity is determined by using a linear gap aligner with 1, -3, -2 for match, mismatch, and gap scores. If the resulting alignment score normalized by the aligned haplotype length is at least 0.97 then the haplotypes are treated as duplicates. Within each duplicate haplotype group, one member is chosen as representative based on having the highest supporting read count from contig assembly, or longest contig size when supporting read count is tied.

### Genotyping
SVs are genotyped in the context of the overlapping haplotype pools created for the purpose of duplicate haplotype identification. Each breakend of each SV allele is evaluated by aligning segments of locally mapped reads to the corresponding segment of the SV haplotype assembly, in addition to aligning these reads to the reference sequence and other SV alleles in the overlapping haplotype pool. The segments of the read, reference and SV haplotypes selected for this purpose on each candidate SV breakend are extended up to 500 bases from the non-anchored edge of the breakend’s homology range. Alignment is scored using a linear gap aligner with 1, -3, -2 for match, mismatch, and gap scores, with alignment scores normalized by the aligned read segment length.

The alignment scores for each read are next converted into support counts, where each read could be identified as uniquely consistent with either the haplotype of the SV allele in question, the reference haplotype or a candidate overlapping SV haplotype. Candidate overlapping haplotypes are only identified for local indel-like SV candidates. The overlapping case can be from a second SV haplotype assembled from reads in the given sample. If a second haplotype wasn’t assembled for the sample, then overlapping SV haplotypes from other samples (among those remaining after merging duplicate haplotypes) are compared to find the overlapping SV haplotype with the highest support from all sample reads. If such an overlapping case is found it is added as a second ‘guest’ candidate haplotype.

Each read’s support for the SV haplotype is evaluated separately for each breakend of a given breakpoint. A read is counted as supporting an SV allele if its alignment score to the SV haplotype is better than the reference haplotype for any breakend. A read can only support the reference haplotype if the reference haplotype alignment has a higher score on all evaluated breakends. This arrangement helps to counteract various forms of reference bias in the alignment scores. Reads supporting an overlapping SV haplotype are internally recorded for the overlapping haplotype but converted into reference allele support for both quality score calculations and in the final VCF allele count output, per the reporting convention used by most available SV calling tools.

### Quality model
The sawfish quality scores are generated from a simplified model which generates qualities directly from the read support counts for each allele. Ideally the scores of read alignments to each allele haplotype would be used instead of counts, but thus far the count approximation hasn’t been attributed to any substantial fraction of SV genotyping errors.

For each SV allele and each sample, we solve for the posterior probability of diploid genotypes, $G$, given the observed read counts at the SV loci, $D$ as follows

$$
P(G \vert D) \propto P(D \vert G) P(G)
$$

The diploid genotypes are $G = \{\text{ref}, \text{het}, \text{hom}\}$ representing 0, 1 or 2 copies of the SV allele. The genotype prior is

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

where $A = \lbrace\text{ref}, \text{alt}\rbrace$ are the SV alleles in the model representing the reference and SV haplotypes. As previously discussed, support for any overlapping SV haplotypes is counted towards the reference allele for the purpose of the genotype quality model, which is intentionally simplified to represent only one SV haplotype at a time.

Considering the terms in the read likelihood, the allele likelihood $P(d \vert a)$ is set from the read allele support counts using a single erroneous read support count probability $e = 1 \times 10^{-5}$ for all cases. The allele probabilities $P(a \vert G)$ are the simple allele fractions (0, 0.5, 1) associated with each genotype.

### Evaluating read coverage support for larger SV calls
After the above described read-based genotyping steps, all deletion and duplication candidates of at least 50kb are additionally evaluated for an SV depth signature consistent with the SV type. If the supporting SV depth signature is not found, the corresponding SV breakpoint is still reported in the VCF output, but as a set of breakend (BND) records, rather than as a deletion or duplication. In this way, breakpoints comprising larger multi-breakpoint complex SVs can still be reported in the VCF output without compromising the precision of simpler copy-number changing deletion and tandem duplication events.

The criteria for a consistent depth signature in each sample is based on the average depth of the interior region of the SV is compared to the average depth in the 6kb regions to the left and right flank. The ratio of the interior to flank depth must be no more than 0.8 for a deletion and no less than 1.2 for a duplication. In a multi-sample context, a consistent depth signature is required in at least half of samples with a non-reference genotype to retain the SV type as a deletion or duplication. While this approach applies a very simple heuristic to assess SV support for each large SV, this has been found to be effective in filtering typical complex SV breakpoints from true large deletion and duplication signatures.
