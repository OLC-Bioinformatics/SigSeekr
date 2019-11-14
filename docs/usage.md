# Quickstart

If you want to run SigSeekr right away upon [installing](installation.md) it, you can do so with a toy dataset.

This dataset is hosted on figshare - to get it, run the following command:

- `wget https://ndownloader.figshare.com/files/9885379 && tar xf 9885379`

You should now have a folder called `example-data` in your present working directory. To run SigSeekr, enter the following command:

- `sigseekr.py -i example-data/inclusion/ -e example-data/exclusion/ -o sigseekr_output -pcr -p3` 

The directory specified with the `-o` flag can be anything - it's the name of a directory where the output files will be created.
Upon entering the command, you should see output that is something like this:

```bash
2019-11-14 10:59:09 Creating inclusion kmer set...
2019-11-14 10:59:25 Creating exclusion kmer set...
2019-11-14 10:59:52 Subtracting exclusion kmers from inclusion kmers with cutoff 1...
2019-11-14 10:59:54 Found kmers unique to inclusion...
2019-11-14 10:59:57 Generating contiguous sequences from inclusion kmers...
2019-11-14 11:00:41 Generating PCR info...
2019-11-14 11:00:43 Finding amplicons of size 200...
2019-11-14 11:01:26 Running Primer3 on potential amplicons...
2019-11-14 11:01:31 SigSeekr run complete!

```

The `sigseekr_output` folder should have five files in it: 

- `amplicons.csv` : list of primers predicted by primer3 and the sizes of their products
- `inclusion_kmers.fasta`: lists all the kmers that are unique to the inclusion set
- `sigseekr_log.txt`: logfile of captured STDOUT and STDERR strings 
- `confirmed_amplicons_200.fasta`: FASTA-formatted file of amplicons present in all inclusion genomes (200 refers to 
the amplicon size specified in the arguments. Default is 200, if multiple sizes are desired, multiple versions of
this file will be created.)
- `potential_pcr_200.fasta`: all potential amplicons based on user-specified amplicon length. This file will be further 
refined by the filtering of amplicon sequences present in the exclusion genomes
- `sigseekr_result.fasta`: regions that unique kmers span


The `sigseekr_result.fasta` created by running SigSeekr on this toy dataset will have one unique region. 

```bash
>contig1_sequence1
AACAGGCGACAGGCAGCATCACTAGCTACTA
```

### Detailed Usage

Detailed usage options can be found by typing `sigseekr.py --help`, which will give the following output. 
Further details on each option can be found below.

```
usage: sigseekr.py [-h] -i INCLUSION -e EXCLUSION -o OUTPUT_FOLDER
                   [-s KMER_SIZE] [-t THREADS] [-pcr] [-k]
                   [-p PLASMID_FILTERING] [-l] [-p3]
                   [-a AMPLICON_SIZE [AMPLICON_SIZE ...]]
                   [-m MAX_POTENTIAL_AMPLICONS]

optional arguments:
  -h, --help            show this help message and exit
  -i INCLUSION, --inclusion INCLUSION
                        Path to folder containing genome(s) you want signature sequences for. Genomes can be in FASTA 
                        or FASTQ format. FASTA-formatted files should be uncompressed, FASTQ-formatted files can be 
                        gzip-compressed or uncompressed.
  -e EXCLUSION, --exclusion EXCLUSION
                        Path to folder containing exclusion genome(s) - those you do not want signature sequences for. 
                        Genomes can be in FASTA or FASTQ format. FASTA-formatted files should be uncompressed, 
                        FASTQ-formatted files can be gzip-compressed or uncompressed.
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Path to folder where you want to store output files. Folder will be created if it does not 
                        exist.
  -s KMER_SIZE, --kmer_size KMER_SIZE
                        Kmer size used to search for sequences unique to inclusion. Default 31. 
                        No idea how changing this affects results. TO BE INVESTIGATED.
  -t THREADS, --threads THREADS
                        Number of threads to run analysis on. Defaults to number of cores on your machine.
  -pcr, --pcr           Enable to filter out inclusion kmers that have close relatives in exclusion kmers.
  -k, --keep_tmpfiles   If enabled, will not clean up a bunch of (fairly) useless files at the end of a run.
  -p PLASMID_FILTERING, --plasmid_filtering PLASMID_FILTERING
                        To ensure unique sequences are not plasmid-borne, a FASTA-formatted database can be provided 
                        with this argument. Any unique kmers that are in the plasmid database will be filtered out.
  -l, --low_memory      Activate this flag to cause plasmid filtering to use substantially less RAM (and go faster), 
                        at the cost of some sensitivity.
  -p3, --primer3        If enabled, will run primer3 on your potential amplicons and generate a list of primers and the 
                        sizes of their products. This output will be found in a file called amplicons.csv in the output 
                        directory specified.
  -a AMPLICON_SIZE [AMPLICON_SIZE ...], --amplicon_size AMPLICON_SIZE [AMPLICON_SIZE ...]
                        Desired size for PCR amplicons. Default 200. If you want to find more than one amplicon size, 
                        enter multiple, separated by spaces.
  -m MAX_POTENTIAL_AMPLICONS, --max_potential_amplicons MAX_POTENTIAL_AMPLICONS
                        If inclusion sequences are very different from exclusion sequences, amplicon generation can take 
                        forever. Set the number of potential amplicons with this option (default 200)

```
