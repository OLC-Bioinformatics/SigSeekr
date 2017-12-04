# System Requirements

SigSeekr has been tested on Debian-based systems (in particular, Ubuntu and Mint), and should have no issue on other Linux-based distributions.
Though not tested, SigSeekr should also work on MacOSX. Windows is not supported at this time.

SigSeekr should be able to run on machines with as little as 8GB of RAM, provided that the `--low_memory` flag is enabled. It is also recommended that a decent amount of disk space is free (100GB for large runs), as the temporary files created in the kmer counting steps in the pipeline can use quite a bit of disk space.

Any number of threads is usable with SigSeekr, with more generally being better.

# Installing External Dependencies

To run SigSeekr, you'll need to install a number of external programs the pipeline uses, and add them to your $PATH.
The programs SigSeekr needs installed are:

- [BBTools >= 37.23](https://jgi.doe.gov/data-and-tools/bbtools/)
- [KMC >= 3.0](http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc&subpage=download)
- [Bedtools >= 2.25.0](https://github.com/arq5x/bedtools2/releases/)
- [Samtools >= 1.6](http://www.htslib.org/download/)

Instructions on how to add a program to your $PATH can be found [here](https://stackoverflow.com/questions/14637979/how-to-permanently-set-path-on-linux-unix).

# Installation via Pip

The easiest way to get SigSeekr up and running is by installing via pip. It's recommended that you create a virtual environment first, and then install.
To create and activate a virtualenv (with python3) something like the following set of commands should work:

- `mkdir SigSeekr`
- `virtualenv -p /usr/bin/python3 SigSeekr`
- `source SigSeekr/bin/activate`

More instructions on virtualenv creation and why virtual environments are wonderful can be found [here](https://realpython.com/blog/python/python-virtual-environments-a-primer/). Once inside the virtual environment, all you should need to do is run pip install:

- `pip install sigseekr`

This command should install any necessary python package dependencies in your virtual environment, and make the SigSeekr script accessible from your terminal. You should now be able to type `sigseekr.py -h` into your terminal and have the help menu for SigSeekr come up.

# Installation from Source

You can also download SigSeekr from GitHub. Releases that are not stable may be pushed to GitHub, so be careful. 
To clone the GitHub repository, type the following:

- `git clone https://github.com/lowandrew/SigSeekr.git`

This should create a folder called `SigSeekr` in your present working directory. You'll need to install any python packages that SigSeekr needs. Within the SigSeekr folder, you should find a file called `requirements.txt`. You can use pip to install all the dependencies that SigSeekr needs with `pip install -r requirements.txt`. If you aren't within a virtual environment, you'll probably need to add a `sudo` before the pip install command.

The `sigseekr.py` script resides in a directory called `sigseekr`. You can either run SigSeekr from inside that directory, or add that directory to your $PATH to have SigSeekr accessible from anywhere.


