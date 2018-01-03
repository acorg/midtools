# Multiple Infection Detection

The various Python scripts in this directory are as follows

## agree-disagree.py

Finds reads that agree with one another at a set of significant genome
locations, makes a graph with reads as nodes and edges between reads that
agree with one another, then finds the connected components of that graph.
The idea / hope was to be able to split the reads into mutually-supporting
groups.

The end result, so far, is that while the reads can be partitioned in this
way the large partitions do not overlap.

## consistency-basic.py

Produces a simple plot showing Adjusted Rand Index (ARI) and Adjusted
Mutual Information (AMI) numbers at a set of significant genome locations.

## consistency-heatmap.py

Consider pairs from a set of significant genome locations. For each pair,
get all reads that cover both locations in the pair. Each location splits
the reads into 4 categories (ACGT). Compare the two splits using ARI and
Normalized Information Distance (NID). These scores are plotted in a
heatmap.

## coverage-and-significant-locations.py

Produces two simple plots. Read coverage level across the genome and the
position of the significant locations.

## significant-base-frequencies.py

For each position in a set of significant genome locations, get the base
frequencies at the position and sort them. Plot a vertical bar for each
location. This looks like a stacked bar chart. It shows sorted base
frequencies though, not the actual bases. The idea is to be able to look to
see if there are many locations that have more than one base present in
significant numbers.

## create-reads.py

Create reads, sampled from a given genome, optionally aligned and mutated.

## extract-reads.py

Extract the reads from a consensus. Removes the gaps and drops the consensus.

## mutate-reads.py

Read a set of reads and mutate them with a given rate.
