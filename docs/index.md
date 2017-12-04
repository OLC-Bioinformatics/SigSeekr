# What is SigSeekr?

SigSeekr is a python pipeline for finding regions unique to one (or multiple) genomes when compared against an 
exclusion set of genomes. These regions can be of use as probes (either bioinformatically, or in a lab-based PCR approach)
to quickly identify whether or not a new genome belongs with the genomes from the inclusion group used to generate the probes.

## How does SigSeekr Work?

SigSeekr uses a kmer-based approach for identifying regions that are unique to an inclusion group. It will create a list of all kmers that are (with size k=31)
that are common across all inclusion genomes specified, and then generates a list of all kmers that are found in the specified exclusion genomes. Any inclusion kmers that
have exact matches to exclusion kmers are then removed, leading to a list of kmers that have at least 1 nucleotide of difference from any exclusion kmer. In the event that
no kmers unique to inclusion sequences are found the process will be repeated, but this time ignoring exclusion kmers with only one occurrence. (This cycle will continue if no
unique kmers are found at one occurrence, going to two, then three, etc.)

The kmers found in the above steps will be suitable for bioinformatic purposes, but may not be sufficiently different from exclusion kmers to be used in a PCR reaction. In order 
to find PCR-appropriate kmers, and additional filtering step can be performed afterwards to only take kmers that are at least 3 SNPs different from any exclusion kmers.
