# Terminology

* Reference - a genome that reads are aligned to.
* Alignment (produced by BWA, Bowtie, BLAST, etc), often as a SAM/BAM file.
* Sites - locations on the genome.
* Significant sites - sites with significant variation in the alignment.

# Be careful

* Some aligners (e.g., `bwa mem`) only align each read to its best-matching
  reference.

What assumptions are you making?

* In simulations, you're making independent reads from a single genome. So
  SNPs are completely uncorrelated. Won't be like that in a real infection.

# Single infections

* What does a single infection with zero mutation rate look like?
* What does a single infection with a non-zero mutation rate look like?
* What does a double infection (with two close viral genomes) mapped
  against one of the genomes look like?

# Multiple Infection Detection

* How can you tell when you have a double (or more) infection?
* Must align against something - else we're reduced to assembly.

## What problems are you trying to solve?

1. When there is a multiple infection of two reasonably similar viruses, a
   matching algorithm will align reads from both to both genomes. A
   consensus made from such a matching is likely to be wrong (depends on
   the relative quantity of virions and how the sequencing process may
   alter the balance of reads in its output).

1. Should be able to pull a clean genome out of mixed reads when only one
   reference is available.  This is a bit like removing noise (e.g., due to
   damage).

1. Given two references, should be able to extract two consensuses.

1. Given a mixed infection (of genetically similar viruses) but only one
   reference genome, should be able to extract a consensus based on and
   alignment to the reference but also a second consensus corresponding to
   the unknown genome.

1. Why can't an algorithm like BWA or BLAST do this already? What are we
   adding? The subtraction of reads that don't match the consensus of the
   reads that do match.  It's not the same as throwing away reads that
   don't match the reference - that would be relatively easy, you'd just
   make the matcher more strict.

1. There is some ratio of virions from one virus to virions from the
   other. This is eventually reflected in some ratio of reads from one or
   the other. The two viruses may be very similar (or identical) in parts
   and different in others. For the identical regions it may not be
   possible (and of course it doesn't matter) to tell which virus a read
   came from.

The various Python scripts in this directory are as follows

## base-frequencies.py

Given a file of aligned reads, print the nucleotide frequencies at each
location.

## connected-components.py

Finds reads that agree with one another at a set of significant genome
locations, makes a graph with reads as nodes and edges between reads that
agree with one another, then finds the connected components of that graph.
The idea is to split the reads into mutually-supporting groups.

## consistency-basic.py

Produces a simple plot showing Adjusted Rand Index (ARI) and Adjusted
Mutual Information (AMI) numbers at a set of significant genome locations.

This is not particularly informative, I don't think.

## consistency-heatmap.py

Consider pairs from a set of significant genome locations. For each pair,
get all reads that cover both locations in the pair. Each location splits
the reads into 4 categories (ACGT). Compare the two splits using ARI and
Normalized Information Distance (NID). The scores are plotted in a heatmap.

Also not very informative.

## coverage-and-significant-locations.py

Produces two simple plots. Read coverage level across the genome and the
position of the significant locations.

## create-reads.py

Create reads, sampled from a given genome, optionally aligned and mutated.

## entropy-for-frequencies.py

Given a list of label frequencies, print their entropies.

## multiple-significant-base-frequencies.py

Plot multiple sorted significant genome location nucleotide frequencies for
a set of aligned reads.

This processes the JSON files containing the sorted location values
produced using the `--valuesFile` option when running
`significant-base-frequencies.py`.

## mutate-reads.py

Read a set of reads and mutate them with a given rate.

## significant-base-frequencies.py

For each position in a set of significant genome locations, get the base
frequencies at the position and sort them. Plot a vertical bar for each
location. This looks like a stacked bar chart. It shows sorted base
frequencies, not the actual bases. The idea is to be able to look to see if
there are many locations that have more than one base present in
significant numbers.
