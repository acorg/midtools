# Multiple Infection Detection

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

## extract-alignment-reads.py

Extract the reads aligned to a consensus, removing gaps. Note that this can
now be done using `filter-fasta.py --idLambda 'lambda id: id.strip("-")'`.

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
frequencies though, not the actual bases. The idea is to be able to look to
see if there are many locations that have more than one base present in
significant numbers.
