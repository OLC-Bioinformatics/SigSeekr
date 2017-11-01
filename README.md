# USS-Ppip
Unique Strain-Specific Python Probe Identifying Program

## Requirements
- mash (https://github.com/marbl/mash) installed and present in your $PATH.
- Python3.5
- BLASTn
- Various python packages (see requirements.txt)
- BBTools (>= 37.23 - https://jgi.doe.gov/data-and-tools/bbtools/) installed and on your $PATH.

## Purpose
Allows the user to find probes that will differentiate between closely related genomes.

## Use
Core.py is the main handler for typing and unique sequence searching methods
* Output folder where a csv typing file and Unique sequences folder will be placed `(-o, --output)`
* Input folder containing genomes you want to differentiate from (in *.fasta format) `(-i, --input)`
* Target genome that you want to make probes for. If this is a folder, all genomes in it will be concatenated into a single genome. `(-t, --target)`

#### Optional Parameters
* Evalue start (default 1e-40) `(-e, --evalue)`
* Evalue stop (default 1e-90) `(-s, --estop)`
* Mash cutoff to eliminate non-target genomes that are extremely similar to one another. (default 0.0002 - higher means more genomes eliminated) `(-c, --mash_cutoff)`
* Number of threads to run mash analyses with. Default is number of CPUs on your system. `(-n, --num_threads)`

#### Example
`Core.py -o USS_Ppip/ -i input_genomes/ -t Ecoli_EDL933.fa`
