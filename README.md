[![Build status](https://travis-ci.org/lowandrew/SigSeekr.svg?master)](https://travis-ci.org/lowandrew)


# SigSeekr

SigSeekr is now new and improved - it uses kmers to find sequences in a set of inclusion genomes that are not present in an exclusion set.

### Installation

Just clone this repository - you'll need the `SigSeekr.py` script.

### External Dependencies

To run SigSeekr, you will need to have the following external programs installed and present on your $PATH:
- BBTools >= 37.23 (https://jgi.doe.gov/data-and-tools/bbtools/)
- KMC >= 3.0 (http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc&subpage=download)
- Python >= 3.5
- bedtools >= 2.25.0 (https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz)
- samtools >= 1.6 (http://www.htslib.org/download/)
 
### Python Package Dependencies

Included in requirements.txt - to install, use pip: `pip install -r requirements.txt`

### Usage

SigSeekr needs a folder containing genomes you want signature sequences (`--inclusion`) for and a folder
containing genomes you do _not_ want the signature sequences matching to (`--exclusion`)

The output folder you specify will end up up containing the following files:
- inclusion_kmers.fasta: A fasta-formatted file containing kmers that are present in inclusion genomes, but not exclusion.
- sigseekr_result.fasta: The kmers in inclusion_kmers mapped back to one of the FASTA-formatted inclusion
genomes to find larger, contiguous unique regions. If you provided only FASTQ files in the
inclusion folder, this file will not be generated.

If you specified the `-pcr` option to try to find kmers suitable for use as PCR primers, the following
 files will also have been created:
 - pcr_kmers.fasta: FASTA-formatted file containing inclusion kmers that have no close matches to 
 any exclusion kmers.
 - amplicons.csv: A comma-separated file showing the potential amplicon size for each possible combination of
 primer kmers. If you provided only FASTQ files in the inclusion folder, this file will not be generated.


```
usage: SigSeekr.py [-h] -i INCLUSION -e EXCLUSION -o OUTPUT_FOLDER
                   [-t THREADS] [-pcr] [-k]

optional arguments:
  -h, --help            show this help message and exit
  -i INCLUSION, --inclusion INCLUSION
                        Path to folder containing genome(s) you want signature
                        sequences for. Genomes can be in FASTA or FASTQ
                        format. FASTA-formatted files should be uncompressed,
                        FASTQ-formatted files can be gzip-compressed or
                        uncompressed.
  -e EXCLUSION, --exclusion EXCLUSION
                        Path to folder containing exclusion genome(s) - those
                        you do not want signature sequences for. Genomes can
                        be in FASTA or FASTQ format. FASTA-formatted files
                        should be uncompressed, FASTQ-formatted files can be
                        gzip-compressed or uncompressed.
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Path to folder where you want to store output files.
                        Folder will be created if it does not exist.
  -t THREADS, --threads THREADS
                        Number of threads to run analysis on. Defaults to
                        number of cores on your machine.
  -pcr, --pcr           Enable to filter out inclusion kmers that have close
                        relatives in exclusion kmers.
  -k, --keep_tmpfiles   If enabled, will not clean up a bunch of (fairly)
                        useless files at the end of a run.

```
