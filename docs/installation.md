# System Requirements

SigSeekr has been tested on Debian-based systems (in particular, Ubuntu and Mint), and should have no issue on other Linux-based distributions.
Though not tested, SigSeekr should also work on MacOSX. Windows is not supported at this time.

SigSeekr should be able to run on machines with as little as 8GB of RAM, provided that the `--low_memory` flag is enabled. It is also recommended that a decent amount of disk space is free (100GB for large runs), as the temporary files created in the kmer counting steps in the pipeline can use quite a bit of disk space.

Any number of threads is usable with SigSeekr, with more generally being better.

# Installation via conda

The easiest way to get SigSeekr up and running is by installing via conda. It's recommended that you create a conda environment first, and then install.

If you need to download and install miniconda:

```bash
$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
$ bash miniconda.sh -b -p $HOME/miniconda
$ export PATH="$HOME/miniconda/bin:$PATH"
$ conda install python=3.6.8
$ conda update -q conda
``` 

If you need to add the required conda-forge and bioconda channels:

```bash  
$ conda config --add channels conda-forge
$ conda config --add channels bioconda
```

__Recommended: create a conda environment:__

```bash
$ conda create -n sigseekr python=3.6.8
$ conda activate sigseekr
```

Install the SigSeekr package from the olcbioinformatics conda channel

```bash
$ conda install -y -c olcbioinformatics sigseekr=0.2.3=py_0
```

This command should install all dependencies, and make the SigSeekr script accessible from your terminal. 
You should now be able to type `sigseekr.py -h` into your terminal and have the help menu for SigSeekr come up.

You can now [run](usage.md) SigSeekr.
