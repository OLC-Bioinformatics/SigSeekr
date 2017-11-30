#!/usr/bin/env python

import multiprocessing
import subprocess
import argparse
import textwrap
import shutil
import glob
import time
import os
import re
import pysam
from Bio import SeqIO
from biotools import kmc
from biotools import bbtools
from accessoryFunctions.accessoryFunctions import printtime


class PcrInfo:
    """
    Class to keep track of things when determining PCR amplicon sizes.
    """
    def __init__(self, sequence, start, end, contig):
        self.seq = sequence
        self.start_position = start
        self.end_position = end
        self.contig_id = contig


def find_paired_reads(fastq_directory, forward_id='_R1', reverse_id='_R2'):
    """
    Looks at a directory to try to find paired fastq files. Should be able to find anything fastq.
    :param fastq_directory: Complete path to directory containing fastq files.
    :param forward_id: Identifier for forward reads. Default R1.
    :param reverse_id: Identifier for reverse reads. Default R2.
    :return: List containing pairs of fastq files, in format [[forward_1, reverse_1], [forward_2, reverse_2]], etc.
    """
    pair_list = list()
    fastq_files = glob.glob(os.path.join(fastq_directory, '*.f*q*'))
    for name in fastq_files:
        if forward_id in name and os.path.isfile(name.replace(forward_id, reverse_id)):
            pair_list.append([name, name.replace(forward_id, reverse_id)])
    return pair_list


def find_unpaired_reads(fastq_directory, forward_id='_R1', reverse_id='_R2'):
    """
    Looks at a directory to try to find unpaired fastq files.
    :param fastq_directory: Complete path to directory containing fastq files.
    :param forward_id: Identifier for forward reads. Default R1.
    :param reverse_id: Identifier for reverse reads. Default R2
    :return: List of unpaired fastq files.
    """
    unpaired_list = list()
    fastq_files = glob.glob(os.path.join(fastq_directory, '*.f*q*'))
    for name in fastq_files:
        if forward_id in name and not os.path.isfile(name.replace(forward_id, reverse_id)):
            unpaired_list.append(name)
        elif forward_id not in name and reverse_id not in name:
            unpaired_list.append(name)
        elif reverse_id in name and not os.path.isfile(name.replace(reverse_id, forward_id)):
            unpaired_list.append(name)
    return unpaired_list


def make_inclusion_kmerdb(inclusion_folder, output_db, forward_id='_R1', reverse_id='_R2', tmpdir='tmpinclusion',
                          maxmem='12', threads='2'):
    """
    Given an folder containing some genomes, finds kmers that are common to all genomes, and writes them to output_db.
    Genomes can be in fasta (uncompressed only? check this) or fastq (gzip compressed or uncompressed) formats.
    Kmers found are 31-mers.
    :param inclusion_folder: Path to folder containing your genomes.
    :param output_db: Base name for the kmc database that will be created.
    :param forward_id: Forward read identifier.
    :param reverse_id: Reverse read identifier.
    :param tmpdir: Directory where temporary databases and whatnot will be stored. Deleted upon method completion.
    :param maxmem: Maximum amount of memory to use when kmerizing, in GB.
    :param threads: Number of threads to use. Counterintuitively, should be a string.
    """
    # Make the tmpdir, if it doesn't exist already.
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)
    # Get lists of everything - fasta, paired fastq, unpaired fastq.
    fastas = glob.glob(os.path.join(inclusion_folder, '*.f*a'))
    paired_fastqs = find_paired_reads(inclusion_folder, forward_id=forward_id, reverse_id=reverse_id)
    unpaired_fastqs = find_unpaired_reads(inclusion_folder, forward_id=forward_id, reverse_id=reverse_id)
    # Make a database for each item in each list, and place it into the tmpdir.
    i = 1
    for fasta in fastas:
        kmc.kmc(fasta, os.path.join(tmpdir, 'database{}'.format(str(i))), fm='', m=maxmem, t=threads,
                tmpdir=os.path.join(tmpdir, str(time.time()).split('.')[0]))
        i += 1
    for pair in paired_fastqs:
        kmc.kmc(forward_in=pair[0], reverse_in=pair[1], database_name=os.path.join(tmpdir, 'database{}'.format(str(i))),
                min_occurrences=2,  # For fastqs, make min_occurrence two to hopefully filter out sequencing errors.
                m=maxmem, t=threads, tmpdir=os.path.join(tmpdir, str(time.time()).split('.')[0]))
        i += 1
    for fastq in unpaired_fastqs:
        kmc.kmc(forward_in=fastq, database_name=os.path.join(tmpdir, 'database{}'.format(str(i))),
                min_occurrences=2,  # For fastqs, make min_occurrence two to hopefully filter out sequencing errors.
                m=maxmem, t=threads, tmpdir=os.path.join(tmpdir, str(time.time()).split('.')[0]))
        i += 1
    # Create a command file to allow kmc to get an intersection of all the inclusion databases created and write to our
    # final inclusion database.
    with open(os.path.join(tmpdir, 'command_file'), 'w') as f:
        f.write('INPUT:\n')
        for j in range(i - 1):
            f.write('set{} = {}\n'.format(str(j + 1), os.path.join(tmpdir, 'database{}'.format(str(j + 1)))))
        f.write('OUTPUT:\n{} = '.format(output_db))
        for j in range(i - 1):
            if j < (i - 2):
                f.write('set{}*'.format(str(j + 1)))
            else:
                f.write('set{}\n'.format(str(j + 1)))
    cmd = 'kmc_tools complex {}'.format(os.path.join(tmpdir, 'command_file'))
    with open(os.path.join(tmpdir, 'asdf.txt'), 'w') as f:
        subprocess.call(cmd, shell=True, stderr=f, stdout=f)
    shutil.rmtree(tmpdir)


def make_exclusion_kmerdb(exclusion_folder, output_db, forward_id='_R1', reverse_id='_R2', tmpdir='tmpexclusion',
                          maxmem='12', threads='2'):
    """
    Given an folder containing some genomes, finds all kmers that are present in genomes, and writes them to output_db.
    Genomes can be in fasta (uncompressed only? check this) or fastq (gzip compressed or uncompressed) formats.
    Kmers found are 31-mers.
    :param exclusion_folder: Path to folder containing your genomes.
    :param output_db: Base name for the kmc database that will be created.
    :param forward_id: Forward read identifier.
    :param reverse_id: Reverse read identifier.
    :param tmpdir: Directory where temporary databases and whatnot will be stored. Deleted upon method completion.
    :param maxmem: Maximum amount of memory to use when kmerizing, in GB.
    :param threads: Number of threads to use. Counterintuitively, should be a string.
    """
    # Make the tmpdir, if it doesn't exist already.
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)
    # Get lists of everything - fasta, paired fastq, unpaired fastq.
    fastas = glob.glob(os.path.join(exclusion_folder, '*.f*a'))
    paired_fastqs = find_paired_reads(exclusion_folder, forward_id=forward_id, reverse_id=reverse_id)
    unpaired_fastqs = find_unpaired_reads(exclusion_folder, forward_id=forward_id, reverse_id=reverse_id)
    # Make a database for each item in each list, and place it into the tmpdir.
    i = 1
    for fasta in fastas:
        kmc.kmc(fasta, os.path.join(tmpdir, 'database{}'.format(str(i))), fm='', m=maxmem,
                tmpdir=os.path.join(tmpdir, str(time.time()).split('.')[0]))
        i += 1
    for pair in paired_fastqs:
        kmc.kmc(forward_in=pair[0], reverse_in=pair[1], database_name=os.path.join(tmpdir, 'database{}'.format(str(i))),
                min_occurrences=2, m=maxmem, t=threads,
                tmpdir=os.path.join(tmpdir, str(time.time()).split('.')[0]))
        i += 1
    for fastq in unpaired_fastqs:
        kmc.kmc(forward_in=fastq, database_name=os.path.join(tmpdir, 'database{}'.format(str(i))),
                min_occurrences=2, m=maxmem, t=threads,
                tmpdir=os.path.join(tmpdir, str(time.time()).split('.')[0]))
        i += 1
    # Create a command file to allow kmc to do a union of all the databases you've created and write them to our final
    # exclusion db.
    with open(os.path.join(tmpdir, 'command_file'), 'w') as f:
        f.write('INPUT:\n')
        for j in range(i - 1):
            f.write('set{} = {}\n'.format(str(j + 1), os.path.join(tmpdir, 'database{}'.format(str(j + 1)))))
        f.write('OUTPUT:\n{} = '.format(output_db))
        for j in range(i - 1):
            if j < (i - 2):
                f.write('set{}+'.format(str(j + 1)))
            else:
                f.write('set{}\n'.format(str(j + 1)))
    cmd = 'kmc_tools complex {}'.format(os.path.join(tmpdir, 'command_file'))
    with open(os.path.join(tmpdir, 'asdf.txt'), 'w') as f:
        subprocess.call(cmd, shell=True, stderr=f, stdout=f)
    shutil.rmtree(tmpdir)


def kmers_to_fasta(kmer_file, output_fasta):
    """
    Given a kmer dump created by using kmc.dump on a kmc database, will transform into a fasta-formatted file.
    :param kmer_file: Path to kmer file.
    :param output_fasta: Path to output file.
    """
    with open(kmer_file) as infile:
        lines = infile.readlines()
    with open(output_fasta, 'w') as outfile:
        i = 1
        for line in lines:
            sequence = line.split()[0]  # Sequence is the first thing in the split
            outfile.write('>kmer{}\n'.format(str(i)))
            outfile.write(sequence + '\n')
            i += 1


def remove_n(input_fasta, output_fasta):
    """
    Given a fasta-formatted file with stretches of Ns in it, will give a fasta-formatted file as output that has
    the original fasta file split into contigs, with a split on each string of Ns.
    :param input_fasta: Path to input fasta.
    :param output_fasta: Path to output fasta. Should NOT be the same as input fasta.
    """
    contigs = SeqIO.parse(input_fasta, 'fasta')
    j = 1
    for contig in contigs:
        sequence = str(contig.seq)
        uniques = re.split('N+', sequence)
        with open(output_fasta, 'a+') as outfile:
            i = 1
            for unique in uniques:
                if unique != '':
                    outfile.write('>contig{}_sequence{}\n'.format(str(j), str(i)))
                    unique = textwrap.fill(unique)
                    outfile.write(unique + '\n')
                    i += 1
            j += 1


def replace_by_index(stretch, seq):
    """
    Given a start and end point in a string (in format 'start:end') and a sequence, will replace characters within
    that stretch with the letter N.
    :param stretch: Start and end index to replace (in format 'start:end')
    :param seq: Sequence to change.
    :return: Sequence modified to have Ns where specified by stretch.
    """
    stretch = stretch.split(':')
    start = int(stretch[0])
    end = int(stretch[1])
    seq = seq[:start] + 'N'*(end-start) + seq[end:]
    return seq


def mask_fasta(input_fasta, output_fasta, bedfile):
    """
    Given a bedfile specifying coverage depths, and an input fasta file corresponding to that bedfile, will create
    a new fasta file with 0-coverage regions replace with Ns.
    :param input_fasta: Path to input fasta file.
    :param output_fasta: Path to output fasta file. Should NOT be the same as input fasta file.
    :param bedfile: Bedfile containing coverage depth info.
    """
    to_mask = dict()
    with open(bedfile) as bed:
        lines = bed.readlines()
    for line in lines:
        line = line.rstrip()
        x = line.split()
        coverage = x[-1]
        end = x[-2]
        start = x[-3]
        name = ' '.join(x[:-3])
        if coverage == '0':
            if name in to_mask:
                to_mask[name].append(start + ':' + end)
            else:
                to_mask[name] = [start + ':' + end]
    fasta_in = SeqIO.parse(input_fasta, 'fasta')
    for contig in fasta_in:
        seq = str(contig.seq)
        if contig.description in to_mask:
            for item in to_mask[contig.description]:
                seq = replace_by_index(item, seq)
            with open(output_fasta, 'a+') as outfile:
                outfile.write('>{}\n'.format(contig.description))
                outfile.write(seq + '\n')


def generate_bedfile(ref_fasta, kmers, output_bedfile, tmpdir='bedgentmp', threads='2'):
    """
    Given a reference FASTA file and a fasta-formatted set of kmers, will generate a coverage bedfile for the reference
    FASTA by mapping the kmers back to the FASTA.
    :param ref_fasta: Path to reference FASTA.
    :param kmers: Path to FASTA-formatted kmer file.
    :param output_bedfile: Path to output bedfile.
    :param tmpdir: Temporary directory to store intermediate files. Will be deleted upon method completion.
    :param threads: Number of threads to use for analysis. Must be a string.
    """
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)
    # First, need to generate a bam file - align the kmers to a reference fasta genome.
    bbtools.bbmap(ref_fasta, kmers, os.path.join(tmpdir, 'out.bam'), threads=threads)
    # Once the bam file is generated, turn it into a sorted bamfile so that bedtools can work with it.
    cmd = 'samtools sort {bamfile} -o {sorted_bamfile}'.format(bamfile=os.path.join(tmpdir, 'out.bam'),
                                                               sorted_bamfile=os.path.join(tmpdir, 'out_sorted.bam'))
    subprocess.call(cmd, shell=True)
    # Use bedtools to get genome coverage, so that we know what to mask.
    cmd = 'bedtools genomecov -ibam {sorted_bamfile} -bga' \
          ' > {output_bed}'.format(sorted_bamfile=os.path.join(tmpdir, 'out_sorted.bam'),
                                   output_bed=output_bedfile)
    subprocess.call(cmd, shell=True)
    shutil.rmtree(tmpdir)


def find_primer_distances(primer_file, reference_file, output_dir, tmpdir='primertmp', threads='2'):
    """
    Given a FASTA-formatted file with kmers that might be acceptable to use as primers, will create a CSV file that
    shows the number of base pairs between each set of primers.
    :param primer_file: Path to FASTA-formatted file containing kmers.
    :param reference_file: Path to FASTA-formatted reference file.
    :param output_dir: Directory where output CSV (called amplicons.csv) will be created.
    :param tmpdir: Temporary directory to store intermediate files. Deleted upon completion.
    :param threads: Number of threads to run analysis on. Must be a string.
    """
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)
    # Step 1: Map the primers back to the reference genome. Only allow perfect matches.
    bbtools.bbmap(reference_file, forward_in=primer_file, out_bam=os.path.join(tmpdir, 'out.bam'), ambig='best',
                  perfectmode='true', threads=threads)
    bam_handle = pysam.AlignmentFile(os.path.join(tmpdir, 'out.bam'), 'rb')
    pcr_info_list = list()
    # Step 2: Parse the bamfile to find the sequence, start position, end position, and contig mapped to for each
    # potential primer.
    for match in bam_handle:
        ref_name = bam_handle.getrname(match.reference_id)
        pcr_info_list.append(PcrInfo(match.query_sequence, match.reference_start,
                                     match.reference_end, ref_name))
    # Now iterate through the list of pcr info to find out amplicon distances.
    with open(os.path.join(output_dir, 'amplicons.csv'), 'w') as f:
        f.write('Sequence1,Sequnce2,Amplicon_Size\n')
    for primer_one in pcr_info_list:
        for primer_two in pcr_info_list:
            # To get amplicon distances need to be on same contig and not be the exact same.
            if primer_one.contig_id == primer_two.contig_id and primer_one.seq != primer_two.seq:
                if int(primer_one.end_position) > int(primer_two.end_position):
                    amplicon_size = int(primer_one.end_position) - int(primer_two.start_position)
                else:
                    amplicon_size = int(primer_two.end_position) - int(primer_one.start_position)
                with open(os.path.join(output_dir, 'amplicons.csv'), 'a+') as f:
                    f.write('{primer_one_seq},{primer_two_seq},{amplicon_size}\n'.format(primer_one_seq=primer_one.seq,
                                                                                         primer_two_seq=primer_two.seq,
                                                                                         amplicon_size=str(amplicon_size)))
    shutil.rmtree(tmpdir)


def main(args):
    start = time.time()
    # Make the necessary inclusion and exclusion kmer sets.
    printtime('Creating inclusion kmer set...', start)
    make_inclusion_kmerdb(args.inclusion, os.path.join(args.output_folder, 'inclusion_db'),
                          tmpdir=os.path.join(args.output_folder, 'inclusiontmp'), threads=str(args.threads))
    printtime('Creating exclusion kmer set...', start)
    make_exclusion_kmerdb(args.exclusion, os.path.join(args.output_folder, 'exclusion_db'),
                          tmpdir=os.path.join(args.output_folder, 'exclusiontmp'), threads=str(args.threads))
    # Now start trying to subtract kmer sets, see how it goes.
    exclusion_cutoff = 1
    while exclusion_cutoff < 10:
        printtime('Subtracting exclusion kmers from inclusion kmers with cutoff {}...'.format(str(exclusion_cutoff)), start)
        kmc.subtract(os.path.join(args.output_folder, 'inclusion_db'), os.path.join(args.output_folder, 'exclusion_db'),
                     os.path.join(args.output_folder, 'unique_to_inclusion_db'), exclude_below=exclusion_cutoff)
        kmc.dump(os.path.join(args.output_folder, 'unique_to_inclusion_db'),
                 os.path.join(args.output_folder, 'unique_kmers.txt'))
        # Now need to check if any kmers are present, and if not, increment the counter to allow a more lax search.
        with open(os.path.join(args.output_folder, 'unique_kmers.txt')) as f:
            lines = f.readlines()
        if lines:
            printtime('Found kmers unique to inclusion...', start)
            break
        exclusion_cutoff += 1
    # Convert our kmers to FASTA format for usage with other programs.
    kmers_to_fasta(os.path.join(args.output_folder, 'unique_kmers.txt'),
                   os.path.join(args.output_folder, 'inclusion_kmers.fasta'))
    # Now that we have kmers that are unique to inclusion sequence, need to map them back to an inclusion genome, if
    # we have one in FASTA format. This will allow us to find unique regions, instead of just kmers.
    if len(glob.glob(os.path.join(args.inclusion, '*.f*a'))) > 0:
        printtime('Generating contiguous sequences from inclusion kmers...', start)
        ref_fasta = glob.glob(os.path.join(args.inclusion, '*.f*a'))[0]
        # Get inclusion kmers into FASTA format.
        generate_bedfile(ref_fasta, os.path.join(args.output_folder, 'inclusion_kmers.fasta'),
                         os.path.join(args.output_folder, 'regions_to_mask.bed'),
                         tmpdir=os.path.join(args.output_folder, 'bedtmp'), threads=str(args.threads))
        mask_fasta(ref_fasta, os.path.join(args.output_folder, 'inclusion_sequence.fasta'),
                   os.path.join(args.output_folder, 'regions_to_mask.bed'))
        remove_n(os.path.join(args.output_folder, 'inclusion_sequence.fasta'),
                 os.path.join(args.output_folder, 'sigseekr_result.fasta'))
    # If we want to find PCR primers, try to filter out any inclusion kmers that are close to exclusion kmers.
    if args.pcr:
        printtime('Generating PCR info...', start)
        # First step is to create fasta of exclusion kmers.
        kmc.dump(os.path.join(args.output_folder, 'exclusion_db'),
                 os.path.join(args.output_folder, 'exclusion_kmers.txt'))
        kmers_to_fasta(os.path.join(args.output_folder, 'exclusion_kmers.txt'),
                       os.path.join(args.output_folder, 'exclusion_kmers.fasta'))
        # Now use bbduk with small kmer size (k=15) to filter out inclusion kmers that have exclusions that are close.
        bbtools.bbduk_filter(reference=os.path.join(args.output_folder, 'exclusion_kmers.fasta'),
                             forward_in=os.path.join(args.output_folder, 'inclusion_kmers.fasta'),
                             forward_out=os.path.join(args.output_folder, 'pcr_kmers.fasta'),
                             k='19', threads=str(args.threads))
        # Next step: Get distances between potential primers by mapping back to a reference (if it exists) and getting
        # distances.
        if len(glob.glob(os.path.join(args.inclusion, '*.f*a'))) > 0:
            ref_fasta = glob.glob(os.path.join(args.inclusion, '*.f*a'))[0]
            find_primer_distances(os.path.join(args.output_folder, 'pcr_kmers.fasta'), ref_fasta, args.output_folder,
                                  tmpdir=os.path.join(args.output_folder, 'pcrtmp'), threads=str(args.threads))
    if not args.keep_tmpfiles:
        printtime('Removing unnecessary output files...', start)
        to_remove = glob.glob(os.path.join(args.output_folder, 'exclusion*'))
        to_remove += glob.glob(os.path.join(args.output_folder, 'unique*'))
        to_remove += glob.glob(os.path.join(args.output_folder, '*kmc*'))
        to_remove += glob.glob(os.path.join(args.output_folder, '*.bed'))
        to_remove += glob.glob(os.path.join(args.output_folder, '*sequence*'))
        for item in to_remove:
            try:
                os.remove(item)
            except FileNotFoundError:  # In case anything was already deleted, don't try to delete it twice.
                pass
    printtime('SigSeekr run complete!', start, '\033[1;32m')


if __name__ == '__main__':
    start = time.time()
    num_cpus = multiprocessing.cpu_count()
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inclusion',
                        type=str,
                        required=True,
                        help='Path to folder containing genome(s) you want signature sequences for.'
                             ' Genomes can be in FASTA or FASTQ format. FASTA-formatted files should be '
                             'uncompressed, FASTQ-formatted files can be gzip-compressed or uncompressed.')
    parser.add_argument('-e', '--exclusion',
                        type=str,
                        required=True,
                        help='Path to folder containing exclusion genome(s) - those you do not want signature'
                             ' sequences for. Genomes can be in FASTA or FASTQ format. FASTA-formatted files should be '
                             'uncompressed, FASTQ-formatted files can be gzip-compressed or uncompressed.')
    parser.add_argument('-o', '--output_folder',
                        type=str,
                        required=True,
                        help='Path to folder where you want to store output files. Folder will be created if it '
                             'does not exist.')
    parser.add_argument('-t', '--threads',
                        type=int,
                        default=num_cpus,
                        help='Number of threads to run analysis on. Defaults to number of cores on your machine.')
    parser.add_argument('-pcr', '--pcr',
                        default=False,
                        action='store_true',
                        help='Enable to filter out inclusion kmers that have close relatives in exclusion kmers.')
    parser.add_argument('-k', '--keep_tmpfiles',
                        default=False,
                        action='store_true',
                        help='If enabled, will not clean up a bunch of (fairly) useless files at the end of a run.')
    args = parser.parse_args()
    if not os.path.isdir(args.output_folder):
        os.makedirs(args.output_folder)
    main(args)
