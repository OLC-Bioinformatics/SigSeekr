# Quickstart

If you want to run SigSeekr right away upon [installing](installation.md) it, you can do so with a toy dataset.

This dataset is hosted on figshare - to get it, run the following command:

- `wget https://ndownloader.figshare.com/files/9885379 && tar xf 9885379`

You should now have a folder called `example-data` in your present working directory. To run SigSeekr, enter the following command:

- `sigseekr.py -i example-data/inclusion/ -e example-data/exclusion/ -o sigseekr_output` 

The directory specified with the `-o` flag can be anything - it's the name of a directory where the output files will be created.
Upon entering the command, you should see output that is something like this:

```bash
 [Elapsed Time: 0.00 seconds] Creating inclusion kmer set... 
 [Elapsed Time: 0.41 seconds] Creating exclusion kmer set... 
 [Elapsed Time: 0.83 seconds] Subtracting exclusion kmers from inclusion kmers with cutoff 1... 
 [Elapsed Time: 0.97 seconds] Found kmers unique to inclusion... 
 [Elapsed Time: 0.97 seconds] Generating contiguous sequences from inclusion kmers... 
 [Elapsed Time: 4.10 seconds] Removing unnecessary output files... 
 [Elapsed Time: 4.10 seconds] SigSeekr run complete! 
```

The `sigseekr_output` folder should have two files in it: `inclusion_kmers.fasta`, which lists all the kmers that are unique to the inclusion set, and `sigseekr_result.fasta`, which contains the regions that unique kmers span. In this case, `sigseekr_result.fasta` should have one unique region. To take a look at it, use the `cat` command:

- `cat sigseekr_output/sigseekr_result.fasta`

The result that should come out of this is:

```bash
>contig1_sequence1
AACAGGCGACAGGCAGCATCACTAGCTACTA
```

### Detailed Usage

Detailed usage options can be found by typing `sigseekr.py --help`, which will give the following output. 
Further details on each option can be found below.

```python
usage: sigseekr.py [-h] -i INCLUSION -e EXCLUSION -o OUTPUT_FOLDER
                   [-t THREADS] [-pcr] [-k] [-p PLASMID_FILTERING] [-l]

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
  -p PLASMID_FILTERING, --plasmid_filtering PLASMID_FILTERING
                        To ensure unique sequences are not plasmid-borne, a
                        FASTA-formatted database can be provided with this
                        argument. Any unique kmers that are in the plasmid
                        database will be filtered out.
  -l, --low_memory      Activate this flag to cause plasmid filtering to use
                        substantially less RAM (and go faster), at the cost of
                        some sensitivity.

```

Additional info:

- `-i, --inclusion`: Not too much to say about this - it's the folder where you'll want to place any genomes that you want signature sequences for. If you place more than one genome here, 
SigSeekr will look for kmers common to all input genomes and develop a signature sequence based on those. These genomes can be FASTA-formatted assemblies (recommended) in which case they must be 
uncompressed, or raw FASTQ reads, in which case they can be either uncompressed or gzip-compressed.
- `-e, --exclusion`: The collection of genomes you _do not_ want your signature sequences to match to. Same file format rules as the inclusion folder. 
- `-o, --output_folder`: The folder where output files will be stored. Created if it doesn't exist. Recommended that you create a new folder for each run, as outputs will be overwritten from previous runs.
- `-t, --threads`: Number of threads to run SigSeekr with. Recommended to leave at the default setting of all cores on your machine, as most programs in the SigSeekr pipeline scale very well with additional threads.
- `-pcr`: Enable to create two additional output files. The first, `pcr_kmers.fasta`, contains inclusion kmers that _should_ be at least 3 SNPs away from any exclusion kmers, making them good potential candidates for PCR primers. The second, `amplicons.csv`, gives a list of all possible pairings of primer candidates, as well as the size of the amplicon that those two primer candidates would create.
- `-k, --keep_tmpfiles`: By default, a number of fairly boring (but sometimes quite large) files are deleted at the end of a run to save on disk space. Specifying this option will keep them around if you want to inspect them more closely. Files that will be kept around with this option specified include the KMC inclusion and exclusion, and unique to inclusion databases (`inclusion_db`, `exclusion_db`, and `unique_to_inclusion_db`), FASTA files of all inclusion kmers (`inclusion_kmers.fasta`), and a bedfile showing coverage of inclusion kmers across one of the inclusion genomes specified (`regions_to_mask.bed`).
- `-p, --plasmid_filtering`: If you're looking for sequences unique to a genome, you probably don't want them on mobile elements that might not be there the next time you look. To help alleviate this potential problem, you can specify the path to a FASTA-formatted database with this option. Any inclusion kmers found in the database will be excluded from further analysis. A relatively extensive plasmid database (~9000 RefSeq plasmids spanning all of Bacteria), can be downloaded with the following command: `https://ndownloader.figshare.com/files/9827323 && tar xf 9827323`.
This will create a folder called `databases` in your present working directory. Within the folder, `plasmid_db.fasta` is the plasmid database.
- `-l, --low_memory`: Using the above-mentioned plasmid database can be memory-intensive. To help alleviate that, add this flag, which will use less memory and go faster, at the cost of some sensitivity.
 

