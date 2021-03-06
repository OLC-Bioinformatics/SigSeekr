#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import dependency_check, SetupLogging
from genewrappers.biotools import bbtools, kmc
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
from io import StringIO
import multiprocessing
import subprocess
import argparse
import textwrap
import primer3
import logging
import shutil
import glob
import time
import os
import re


def write_to_logfile(logfile, out, err, cmd):
    """
    Write stdout and stderr of a system call, as well as the command to a supplied logfile
    :param logfile: Name and path of the logfile
    :param out: stdout from the command
    :param err: stderr from the command
    :param cmd: command run
    """
    with open(logfile, 'a+') as outfile:
        outfile.write('Command used: {}\n\n'.format(cmd))
        outfile.write('STDOUT: {}\n\n'.format(out))
        outfile.write('STDERR: {}\n\n'.format(err))


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


def make_kmerdb(folder, output_db, tmpdir, analysis, forward_id='_R1', reverse_id='_R2',
                maxmem='12', threads='2', logfile=None, k=31, keep=False):
    """
    Given an folder containing some genomes, finds all kmers that are present in genomes, and writes them to output_db.
    Genomes can be in fasta (uncompressed only? check this) or fastq (gzip compressed or uncompressed) formats.
    Kmers found are 31-mers.
    :param folder: Path to folder containing your genomes.
    :param output_db: Base name for the kmc database that will be created.
    :param tmpdir: Directory where temporary databases and whatnot will be stored. Deleted upon method completion.
    :param analysis: type STR: Options are 'inclusion' or 'exclusion'. For an inclusion database, an intersection of
    kmers is calculated, while a union of all kmers is created for the exlusion database
    :param forward_id: Forward read identifier.
    :param reverse_id: Reverse read identifier.
    :param maxmem: Maximum amount of memory to use when kmerizing, in GB.
    :param threads: Number of threads to use. Counterintuitively, should be a string.
    :param logfile: Text file you want commands used, as well as stdout and stderr from called programs, to be logged to
    :param k: Kmer size to use.
    :param keep: Passed argument on whether to keep temporary files
    """
    # Make the tmpdir, if it doesn't exist already.
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)
    # Get lists of everything - fasta, paired fastq, unpaired fastq.
    fastas = glob.glob(os.path.join(folder, '*.f*a'))
    paired_fastqs = find_paired_reads(folder,
                                      forward_id=forward_id,
                                      reverse_id=reverse_id)
    unpaired_fastqs = find_unpaired_reads(folder,
                                          forward_id=forward_id,
                                          reverse_id=reverse_id)
    # Make a database for each item in each list, and place it into the tmpdir.
    i = 1
    for fasta in fastas:
        database_name = os.path.join(tmpdir, 'database{}'.format(str(i)))
        if not os.path.isfile(database_name + '.kmc_pre'):
            out, err, cmd = kmc.kmc(forward_in=fasta,
                                    database_name=database_name,
                                    fm='',
                                    m=maxmem,
                                    t=threads,
                                    tmpdir=os.path.join(tmpdir, str(time.time()).split('.')[0]),
                                    returncmd=True,
                                    k=k)
            if logfile:
                write_to_logfile(logfile=logfile, out=out, err=err, cmd=cmd)
        i += 1
    for pair in paired_fastqs:
        database_name = os.path.join(tmpdir, 'database{}'.format(str(i)))
        if not os.path.isfile(database_name + '.kmc_pre'):

            out, err, cmd = kmc.kmc(forward_in=pair[0],
                                    reverse_in=pair[1],
                                    database_name=database_name,
                                    min_occurrences=2,
                                    m=maxmem,
                                    t=threads,
                                    tmpdir=os.path.join(tmpdir, str(time.time()).split('.')[0]),
                                    returncmd=True,
                                    k=k)
            if logfile:
                write_to_logfile(logfile=logfile, out=out, err=err, cmd=cmd)
        i += 1
    for fastq in unpaired_fastqs:
        database_name = os.path.join(tmpdir, 'database{}'.format(str(i)))
        if not os.path.isfile(database_name + '.kmc_pre'):

            out, err, cmd = kmc.kmc(forward_in=fastq,
                                    database_name=database_name,
                                    min_occurrences=2,
                                    m=maxmem,
                                    t=threads,
                                    tmpdir=os.path.join(tmpdir, str(time.time()).split('.')[0]),
                                    returncmd=True,
                                    k=k)
            if logfile:
                write_to_logfile(logfile=logfile, out=out, err=err, cmd=cmd)
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
                if analysis == 'inclusion':
                    f.write('set{}*sum'.format(str(j + 1)))
                else:
                    f.write('set{}+'.format(str(j + 1)))
            else:
                f.write('set{}\n'.format(str(j + 1)))
    cmd = 'kmc_tools complex {}'.format(os.path.join(tmpdir, 'command_file'))
    if not os.path.isfile(output_db + '.kmc_pre'):
        if logfile:
            with open(logfile, 'a+') as f:
                f.write('Command: {}'.format(cmd))
                subprocess.call(cmd, shell=True, stderr=f, stdout=f)
        else:
            with open(os.devnull, 'w') as f:
                subprocess.call(cmd, shell=True, stderr=f, stdout=f)
    if not keep:
        shutil.rmtree(tmpdir)


def kmers_to_fasta(kmer_file, output_fasta):
    """
    Given a kmer dump created by using kmc.dump on a kmc database, will transform into a fasta-formatted file.
    :param kmer_file: Path to kmer file.
    :param output_fasta: Path to output file.
    """
    if not os.path.isfile(output_fasta):
        with open(kmer_file) as infile:
            lines = infile.readlines()
        with open(output_fasta, 'w') as outfile:
            i = 1
            for line in lines:
                sequence = line.split()[0]  # Sequence is the first thing in the split
                outfile.write('>kmer{}\n'.format(str(i)))
                outfile.write(sequence + '\n')
                i += 1


def remove_n(input_fasta, output_fasta, k=31):
    """
    Given a fasta-formatted file with stretches of Ns in it, will give a fasta-formatted file as output that has
    the original fasta file split into contigs, with a split on each string of Ns.
    :param input_fasta: Path to input fasta.
    :param k: kmer size used to find unique kmers. Unique sequence need to be at least this long.
    :param output_fasta: Path to output fasta. Should NOT be the same as input fasta.
    """
    if not os.path.isfile(output_fasta):
        contigs = SeqIO.parse(input_fasta, 'fasta')
        j = 1
        for contig in contigs:
            sequence = str(contig.seq)
            uniques = re.split('N+', sequence)
            with open(output_fasta, 'a+') as outfile:
                i = 1
                for unique in uniques:
                    if unique != '' and len(unique) >= k:
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
    seq = seq[:start] + 'N' * (end - start) + seq[end:]
    return seq


def mask_fasta(input_fasta, output_fasta, bedfile, k=31):
    """
    Given a bedfile specifying coverage depths, and an input fasta file corresponding to that bedfile, will create
    a new fasta file with 0-coverage regions replace with Ns.
    :param input_fasta: Path to input fasta file.
    :param output_fasta: Path to output fasta file. Should NOT be the same as input fasta file.
    :param bedfile: Bedfile containing coverage depth info.
    :param k: kmer size used when generating kmers unique to inclusion with kmc.
    """
    if not os.path.isfile(output_fasta):
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
            # This needs to be at least k, otherwise you end up with a stretch of length (k * 2) + 1 which matches
            # exclusion sequence perfectly, except for a SNV in the middle.
            if int(coverage) < k:
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


def generate_bedfile(ref_fasta, kmers, output_bedfile, tmpdir='bedgentmp', threads='2', logfile=None, keep=False):
    """
    Given a reference FASTA file and a fasta-formatted set of kmers, will generate a coverage bedfile for the reference
    FASTA by mapping the kmers back to the FASTA.
    :param ref_fasta: Path to reference FASTA.
    :param kmers: Path to FASTA-formatted kmer file.
    :param output_bedfile: Path to output bedfile.
    :param tmpdir: Temporary directory to store intermediate files. Will be deleted upon method completion.
    :param threads: Number of threads to use for analysis. Must be a string.
    :param logfile: Logfile to write stdout and stderr for makeblastdb call to.
    :param keep: Passed argument on whether to keep temporary files
    """
    if not os.path.isfile(output_bedfile):
        if not os.path.isdir(tmpdir):
            os.makedirs(tmpdir)
        # First, need to generate a bam file - align the kmers to a reference fasta genome.
        bbtools.bbmap(ref_fasta, kmers, os.path.join(tmpdir, 'out.bam'), threads=threads, ambig='best',
                      perfectmode='true')
        # Once the bam file is generated, turn it into a sorted bamfile so that bedtools can work with it.
        cmd = 'samtools sort {bamfile} -o {sorted_bamfile}'.format(bamfile=os.path.join(tmpdir, 'out.bam'),
                                                                   sorted_bamfile=os.path.join(tmpdir,
                                                                                               'out_sorted.bam'))
        if logfile:
            with open(logfile, 'a+') as f:
                f.write('Command: {}'.format(cmd))
                subprocess.call(cmd, shell=True, stderr=f, stdout=f)
        else:
            subprocess.call(cmd, shell=True)
        # Use bedtools to get genome coverage, so that we know what to mask.
        cmd = 'bedtools genomecov -ibam {sorted_bamfile} -bga > {output_bed}' \
            .format(sorted_bamfile=os.path.join(tmpdir, 'out_sorted.bam'),
                    output_bed=output_bedfile)
        if logfile:
            with open(logfile, 'a+') as f:
                f.write('Command: {}'.format(cmd))
                subprocess.call(cmd, shell=True, stderr=f, stdout=f)
        else:
            subprocess.call(cmd, shell=True)
        if not keep:
            shutil.rmtree(tmpdir)


def split_sequences_into_amplicons(input_sequence_file, output_amplicon_file, amplicon_length=200,
                                   max_potential_amplicons=200):
    """
    Given an input fasta file, will find potential amplicons of amplicon_length. Uses a sliding-windowy approach,
    so if a contig is 900 bp long and amplicon length is 200, will get 4 amplicons - positions 1-200, 201-400,
    401-600, and 601-800. Last 100 bp will be ignored. If contig is shorter than amplicon_length, nothing happens.
    :param input_sequence_file: Path to input file to be split into amplicons, in FASTA format.
    :param output_amplicon_file: Path to where you'll want output amplicon file to be stored. Will overwrite
    existing files.
    :param amplicon_length: Desired amplicon length. Default 200.
    :param max_potential_amplicons: Maximum number of amplicons to generate. Default is 200.
    """
    if not os.path.isfile(output_amplicon_file):
        with open(output_amplicon_file, 'w') as f:
            seq_id = 1
            for unique_seq in SeqIO.parse(input_sequence_file, 'fasta'):
                if len(unique_seq) >= amplicon_length:
                    sequence = str(unique_seq.seq)
                    outstr = ''
                    # Split sequence into a bunch of chunks of specified amplicon size, write those chunks to file.
                    i = 0
                    while i < len(sequence):
                        if len(sequence[i:amplicon_length + i]) == amplicon_length:
                            outstr += '>sequence' + str(seq_id) + '\n'
                            outstr += sequence[i:amplicon_length + i] + '\n'
                            seq_id += 1
                        i += amplicon_length
                    f.write(outstr)
                if seq_id >= max_potential_amplicons:
                    break


def make_all_exclusion_blast_db(exclusion_folder, combined_exclusion_fasta, logfile=None):
    """
    Given a folder (which we assume contains fasta files) and the name of a desired output file, will concatenate the
    fasta files in the exclusion folder, put them into the combined_exclusion_fasta, and then make a blast database.
    If a logfile is specified, stdout and stderr from the makeblastdb call will be written to it.
    :param exclusion_folder: Path to folder containing fasta files to be concatenated.
    :param combined_exclusion_fasta: Path to desired combined fasta file. Also will be the -db parameter to pass to
    blastn when the time comes for that.
    :param logfile: Logfile to write stdout and stderr for makeblastdb call to.
    """
    if not os.path.isfile(combined_exclusion_fasta):
        # Get all the files copied together with shutil instead of cat to keep things working cross-platform
        with open(combined_exclusion_fasta, 'wb') as wfd:
            for f in glob.glob(os.path.join(exclusion_folder, '*.f*a')):
                with open(f, 'rb') as fd:
                    shutil.copyfileobj(fd, wfd)
    if not os.path.isfile(combined_exclusion_fasta + '.nhr'):
        # With all the files combined, make a blast database, and get stdout and stderr output to logfile (if specified)
        cmd = 'makeblastdb -in {combined_fasta} -dbtype nucl'.format(combined_fasta=combined_exclusion_fasta)
        if logfile:
            with open(logfile, 'a+') as f:
                f.write('Command: {}'.format(cmd))
                subprocess.call(cmd, shell=True, stderr=f, stdout=f)
        else:
            subprocess.call(cmd, shell=True)


def ensure_amplicons_not_in_exclusion(exclusion_blastdb, potential_amplicons, confirmed_amplicons,
                                      max_potential_amplicons=200):
    """
    Given a blast database of sequences we do not want amplicons to match to and a fasta file containing our
    potential amplicons, will blastn potential amplicons to make sure that they don't match too closely to the
    exclusion blastdb. Criteria for this: Top hit length can't be more than 40 base pairs (anything more than
    that might start getting amplified if we're really unlucky) and if more than one hit, can't have any two hits
    within 5000bp of each other, as those could also potentially amplify if we're really unlucky.
    Amplicons confirmed to meet these criteria will get written to confirmed_amplicons, which will overwrite any
    file that was already there.
    :param exclusion_blastdb: Path to exclusion blast database. In this pipeline, should have been created by
    make_all_exclusion_blast_db
    :param potential_amplicons: Path to potential amplicon fasta file. In this pipeline, should have been created by
    split_sequences_into_amplicons
    :param confirmed_amplicons: Path to your desired output confirmed amplicon file. Overwrites file if something
    was already there.
    :param max_potential_amplicons: Maximum number of amplicons to generate. Default is 200.
    """
    if not os.path.isfile(confirmed_amplicons):
        outstr = ''
        sequence_id = 1
        for potential_sequence in SeqIO.parse(potential_amplicons, 'fasta'):
            blastn = NcbiblastnCommandline(db=exclusion_blastdb,
                                           task='blastn',
                                           outfmt=5)
            stdout, stderr = blastn(stdin=str(potential_sequence.seq))
            top_hit_length = 999999  # Start this at ridiculously high value

            # The hit location dict will store the locations of every blast hit to each contig.
            # Each contig is an entry into the dict, with each entry being a list of locations.
            # We'll later try every combination in each list to make sure no two matches are too close together.
            hit_location_dict = dict()
            records = NCBIXML.parse(StringIO(stdout))
            for record in records:
                try:
                    top_hit_length = record.alignments[0].hsps[0].align_length
                except IndexError:  # Should happen if we don't have any hits at all.
                    top_hit_length = 0
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        if alignment.title in hit_location_dict:
                            hit_location_dict[alignment.title].append(hsp.sbjct_start)
                        else:
                            hit_location_dict[alignment.title] = [hsp.sbjct_start]

            # Set up a flag that we'll turn to true if we find any sets of matches that are too close together.
            matches_too_close = False
            for contig in hit_location_dict:
                for i in range(len(hit_location_dict[contig])):
                    for j in range(len(hit_location_dict[contig])):
                        if i != j:
                            # Make sure no two hits within 5000bp of each other.
                            if abs(hit_location_dict[contig][i] - hit_location_dict[contig][j]) < 5000:
                                matches_too_close = True
            # Allow writing to outstr if either we have no hits longer than a roughly two pcr primers (so 40ish bp)
            # Also can't have any two matches to the same contig within 5000bp of each other.
            if top_hit_length < 40 and matches_too_close is False:
                outstr += '>sequence' + str(sequence_id) + '\n'
                outstr += str(potential_sequence.seq) + '\n'
                sequence_id += 1
            if sequence_id > max_potential_amplicons:
                break
        with open(confirmed_amplicons, 'w') as f:
            f.write(outstr)


def confirm_amplicons_in_all_inclusion_genomes(inclusion_fasta_dir, potential_amplicon_file, confirmed_amplicon_file,
                                               tmpdir='tmp', logfile=None, amplicon_size=200, keep=False,
                                               max_potential_amplicons=200):
    """
    Provided with a directory containing fasta files you want your amplicons to match to and a fasta file where each
    entry is a potential amplicon, will ensure that each genome contains a full-length match with at least 99 percent
    identity to each amplicon. Amplicons that do not meet these criteria are filtered out.
    :param inclusion_fasta_dir: Path to directory containing fastas that you want your amplicons to have matches to.
    :param potential_amplicon_file: Potential amplicon fasta file. At this point, should be the file created by
    ensure_amplicons_not_in_exclusion
    :param confirmed_amplicon_file: Path to file where amplicons confirmed to be in all inclusion genomes will be
    written.
    :param tmpdir: Path to directory where blastdbs for each inclusion genome will be created.
    :param logfile: Logfile for stdout and stderr from the makeblastdb commands.
    :param amplicon_size: Desired size to use for amplicon creation. Default is 200.
    :param keep: Passed argument on whether to keep temporary files
    :param max_potential_amplicons: Maximum number of amplicons to generate. Default is 200.
    """
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)
    if not os.path.isfile(confirmed_amplicon_file):
        # Copy each of the inclusion fastas to our temporary folder so we can make them into blast dbs
        inclusion_fastas = glob.glob(os.path.join(inclusion_fasta_dir, '*.f*a'))
        for inclusion_fasta in inclusion_fastas:
            shutil.copy(inclusion_fasta, tmpdir)
        # Now make a blast DB for each. (And write to logfile, if specified!)
        inclusion_fastas = glob.glob(os.path.join(tmpdir, '*.f*a'))
        for inclusion_fasta in inclusion_fastas:
            cmd = 'makeblastdb -in {inclusion_fasta} -dbtype nucl'.format(inclusion_fasta=inclusion_fasta)
            if logfile:
                with open(logfile, 'a+') as f:
                    f.write('Command: {}'.format(cmd))
                    subprocess.call(cmd, shell=True, stderr=f, stdout=f)
            else:
                subprocess.call(cmd, shell=True)
        # Finally, BLAST each of our 'confirmed' amplicons against each database, make sure we get a full length match,
        # which has a percent identity of at least 99 percent in each inclusion genome. If we don't get that, don't
        # write the sequence.
        all_fasta_count = 0
        for potential_sequence in SeqIO.parse(potential_amplicon_file, 'fasta'):
            in_all_fastas = True  # Assume each sequence is in all fastas. Set to False if necessary
            for inclusion_fasta in inclusion_fastas:
                blastn = NcbiblastnCommandline(db=inclusion_fasta,
                                               outfmt=5)
                stdout, stderr = blastn(stdin=str(potential_sequence.seq))
                for record in NCBIXML.parse(StringIO(stdout)):
                    # Attempt to parse top hit info. If, somehow, no hits exist, set our flag to False
                    try:
                        hit_length = record.alignments[0].hsps[0].align_length
                        percent_id = float(record.alignments[0].hsps[0].identities) / float(amplicon_size)
                        if hit_length != amplicon_size or percent_id < 0.99:
                            in_all_fastas = False
                    except IndexError:
                        in_all_fastas = False
            if in_all_fastas:
                all_fasta_count += 1
                with open(confirmed_amplicon_file, 'a+') as f:
                    f.write('>sequence{}\n'.format(all_fasta_count))
                    f.write(str(potential_sequence.seq) + '\n')
            if all_fasta_count >= max_potential_amplicons:
                break
    if not keep:
        shutil.rmtree(tmpdir)


def run_primer3(amplicon_files, output_file):
    """
    Runs primer3 via the primer3-py python binding (which means that the packaging system gets primer3 installed for us,
    which is great). Given a list of fasta files that each contain sequence unique to our inclusion genome, will run
    primer3 on each potential amplicon. Generates a csv file that gives the top 10 (or less if there are fewer) hits
    for each sequence in each file, and tells the user left and right primers, as well as product size.
    :param amplicon_files: List of potential amplicon fasta files.
    :param output_file: Path to output csv file. Will be overwritten if it already exists.
    """
    with open(output_file, 'w') as f:
        f.write('Left_Primer_Sequence,Right_Primer_Sequence,Amplicon_Size\n')
    detailed_amplicons = output_file.replace('.csv', '_detailed.csv')
    best_primers = output_file.replace('.csv', '_best.csv')
    detailed_header = 'SequenceName,Primer,start,len,tm,gc%,any_th,3\'_th,hairpin,seq,penalty\n'
    detailed_body = str()
    primer_dict = dict()
    penalties = set()
    for amplicon_file in amplicon_files:
        primer_dict[amplicon_file] = dict()
        for potential in SeqIO.parse(amplicon_file, 'fasta'):
            primer_dict[amplicon_file][potential.id] = dict()
            # Settings for primer design mostly taken from the defaults on primer3 web implementation found at
            # http://bioinfo.ut.ee/primer3-0.4.0/
            # I'm assuming defaults are generally what most people will want - should do some verification of this.
            primer_info_dict = primer3.bindings.designPrimers({
                'SEQUENCE_ID': potential.id,
                'SEQUENCE_TEMPLATE': str(potential.seq)
            },
                {
                    'PRIMER_OPT_SIZE': 20,
                    'PRIMER_PICK_INTERNAL_OLIGO': 0,
                    'PRIMER_INTERNAL_MAX_SELF_END': 8,
                    'PRIMER_MIN_SIZE': 18,
                    'PRIMER_MAX_SIZE': 25,
                    'PRIMER_OPT_TM': 60.0,
                    'PRIMER_MIN_TM': 57.0,
                    'PRIMER_MAX_TM': 63.0,
                    'PRIMER_MIN_GC': 20.0,
                    'PRIMER_MAX_GC': 80.0,
                    'PRIMER_MAX_POLY_X': 100,
                    'PRIMER_INTERNAL_MAX_POLY_X': 100,
                    'PRIMER_SALT_MONOVALENT': 50.0,
                    'PRIMER_DNA_CONC': 50.0,
                    'PRIMER_MAX_NS_ACCEPTED': 0,
                    'PRIMER_MAX_SELF_ANY': 12,
                    'PRIMER_MAX_SELF_END': 8,
                    'PRIMER_PAIR_MAX_COMPL_ANY': 12,
                    'PRIMER_PAIR_MAX_COMPL_END': 8
                })
            # Now need to parse the giant result dict that primer3 makes
            primer_dict[amplicon_file][potential.id].update(primer_info_dict)
            with open(output_file, 'a+') as f:
                for i in range(primer_info_dict['PRIMER_PAIR_NUM_RETURNED']):
                    product_size = str(primer_info_dict['PRIMER_PAIR_{}_PRODUCT_SIZE'.format(i)])
                    left_primer_seq = primer_info_dict['PRIMER_LEFT_{}_SEQUENCE'.format(i)]
                    right_primer_seq = primer_info_dict['PRIMER_RIGHT_{}_SEQUENCE'.format(i)]
                    f.write('{},{},{}\n'.format(left_primer_seq, right_primer_seq, product_size))
            for i in range(primer_info_dict['PRIMER_PAIR_NUM_RETURNED']):
                # LEFT PRIMER
                left_primer_name = 'PRIMER_LEFT_{num}'.format(num=i)
                detailed_body = populate_detailed_body(body_string=detailed_body,
                                                       primer_name=left_primer_name,
                                                       primer_info_dict=primer_info_dict,
                                                       seq_name=potential.id)
                # RIGHT PRIMER
                right_primer_name = 'PRIMER_RIGHT_{num}'.format(num=i)
                detailed_body = populate_detailed_body(body_string=detailed_body,
                                                       primer_name=right_primer_name,
                                                       primer_info_dict=primer_info_dict)
                # PRIMER PAIR
                detailed_body = populate_body_pair(body_string=detailed_body,
                                                   primer_info_dict=primer_info_dict,
                                                   iterator=i)
                penalties.add(primer_info_dict['PRIMER_PAIR_{}_PENALTY'.format(i)])
    with open(detailed_amplicons, 'w') as details:
        details.write(detailed_header)
        details.write(detailed_body)
    best_penalties = list()
    # Try to find the best five penalties
    try:
        best_penalties = sorted(list(penalties))[:5]
    # If there are fewer than five results returned, add all the penalties
    except IndexError:
        for penalty in sorted(list(penalties)):
            best_penalties.append(penalty)
    # Create a report of the best primer sets
    best_body = str()
    for amplicon_file, seq_name_dict in primer_dict.items():
        for seq_name, primer_info_dict in seq_name_dict.items():
            for i in range(primer_info_dict['PRIMER_PAIR_NUM_RETURNED']):
                if primer_info_dict['PRIMER_PAIR_{}_PENALTY'.format(i)] in best_penalties:
                    # LEFT PRIMER
                    primer_name = 'PRIMER_LEFT_{num}'.format(num=i)
                    best_body = populate_detailed_body(body_string=best_body,
                                                       primer_name=primer_name,
                                                       primer_info_dict=primer_info_dict,
                                                       seq_name=seq_name)
                    # RIGHT PRIMER
                    primer_name = 'PRIMER_RIGHT_{num}'.format(num=i)
                    best_body = populate_detailed_body(body_string=best_body,
                                                       primer_name=primer_name,
                                                       primer_info_dict=primer_info_dict)
                    # PRIMER PAIR
                    best_body = populate_body_pair(body_string=best_body,
                                                   primer_info_dict=primer_info_dict,
                                                   iterator=i)
    with open(best_primers, 'w') as best:
        best.write(detailed_header)
        best.write(best_body)


def populate_detailed_body(body_string, primer_name, primer_info_dict, seq_name=''):
    """
    Populate the primer (forward or reverse) with the appropriate values extracted from the primer dictionary to
    create detailed reports
    :param body_string: type STR: String containing growing primer details
    :param primer_name: type STR: Direction and iterator of primer e.g. PRIMER_LEFT_0
    :param primer_info_dict: type DICT: Dictionary output from primer3
    :param seq_name: type STR: Name of amplicon being tested
    :return: Updated body string
    """
    body_string += '{seqname},{pname},{start},{length},{tm},{gc},{any_th},{end_th},{hairpin},' \
                   '{seq},{pen}\n'\
        .format(seqname=seq_name,
                pname=primer_name,
                start=str(primer_info_dict[primer_name][0] + 1),
                length=str(primer_info_dict[primer_name][1]),
                tm='{:0.2f}'.format(primer_info_dict['{pn}_TM'.format(pn=primer_name)]),
                gc='{:0.2f}'.format(primer_info_dict['{pn}_GC_PERCENT'.format(pn=primer_name)]),
                any_th='{:0.2f}'.format(primer_info_dict['{pn}_SELF_ANY_TH'.format(pn=primer_name)]),
                end_th='{:0.2f}'.format(primer_info_dict['{pn}_SELF_END_TH'.format(pn=primer_name)]),
                hairpin='{:0.2f}'.format(primer_info_dict['{pn}_HAIRPIN_TH'.format(pn=primer_name)]),
                seq=str(primer_info_dict['{pn}_SEQUENCE'.format(pn=primer_name)]),
                pen='{:0.4f}'.format(primer_info_dict['{pn}_PENALTY'.format(pn=primer_name)]))
    return body_string


def populate_body_pair(body_string, primer_info_dict, iterator):
    """
    Populate the primer pair string for the detailed reports
    :param body_string: type STR: String containing growing primer pair details
    :param primer_info_dict: type DICT: Dictionary output from primer3
    :param iterator: type STR: Iterator of number of primer pairs returned (typecast to STR)
    :return: Updated body string
    """
    body_string += \
        ',PRODUCT SIZE:, {ps}, PAIR ANY_TH COMPL:, {compl_any_th}, PAIR 3\'_TH COMPL:, ' \
        '{compl_end_th},PRIMER PAIR PENALTY:, {ppp}\n\n'\
        .format(ps=str(primer_info_dict['PRIMER_PAIR_{}_PRODUCT_SIZE'.format(iterator)] - 1),
                compl_any_th='{:0.2f}'.format(
                    primer_info_dict['PRIMER_PAIR_{}_COMPL_ANY_TH'.format(iterator)]),
                compl_end_th='{:0.2f}'.format(
                    primer_info_dict['PRIMER_PAIR_{}_COMPL_END_TH'.format(iterator)]),
                ppp='{:0.4f}'.format(primer_info_dict['PRIMER_PAIR_{}_PENALTY'.format(iterator)]))
    return body_string


def main(args):
    log = os.path.join(args.output_folder, 'sigseekr_log.txt')
    # Make the necessary inclusion and exclusion kmer sets.
    logging.info('Creating inclusion kmer set...')
    make_kmerdb(folder=args.inclusion,
                output_db=os.path.join(args.output_folder, 'inclusion_db'),
                tmpdir=os.path.join(args.output_folder, 'inclusiontmp'),
                analysis='inclusion',
                threads=str(args.threads),
                maxmem=str(args.max_memory),
                logfile=log,
                k=args.kmer_size,
                keep=args.keep_tmpfiles)
    logging.info('Creating exclusion kmer set...')
    make_kmerdb(folder=args.exclusion,
                output_db=os.path.join(args.output_folder, 'exclusion_db'),
                tmpdir=os.path.join(args.output_folder, 'exclusiontmp'),
                analysis='exclusion',
                threads=str(args.threads),
                maxmem=str(args.max_memory),
                logfile=log,
                k=args.kmer_size,
                keep=args.keep_tmpfiles)
    # Now start trying to subtract kmer sets. If no results are generated with an exclusion cutoff of 1, keep trying
    # with values up to 10. Then give up.
    exclusion_cutoff = 1
    while exclusion_cutoff <= 10:
        logging.info('Subtracting exclusion kmers from inclusion kmers with cutoff {}...'.format(str(exclusion_cutoff)))
        out, err, cmd = kmc.subtract(database_1=os.path.join(args.output_folder, 'inclusion_db'),
                                     database_2=os.path.join(args.output_folder, 'exclusion_db'),
                                     results=os.path.join(args.output_folder, 'unique_to_inclusion_db'),
                                     exclude_below=exclusion_cutoff,
                                     returncmd=True)
        write_to_logfile(logfile=log, out=out, err=err, cmd=cmd)
        # Dump the unique kmers from kmc database format into a text file
        out, err, cmd = kmc.dump(database=os.path.join(args.output_folder, 'unique_to_inclusion_db'),
                                 output=os.path.join(args.output_folder, 'unique_kmers.txt'),
                                 returncmd=True)
        write_to_logfile(logfile=log, out=out, err=err, cmd=cmd)
        # Now need to check if any kmers are present, and if not, increment the counter to allow a more lax search.
        with open(os.path.join(args.output_folder, 'unique_kmers.txt')) as f:
            lines = f.readlines()
        if lines:
            logging.info('Found kmers unique to inclusion...')
            # Convert our kmers to FASTA format for usage with other programs.
            kmers_to_fasta(kmer_file=os.path.join(args.output_folder, 'unique_kmers.txt'),
                           output_fasta=os.path.join(args.output_folder, 'inclusion_kmers.fasta'))
            # Filter out kmers that are plasmid-borne, as necessary.
            if args.plasmid_filtering != 'NA':
                logging.info('Filtering out inclusion kmers that map to plasmids...')
                if args.low_memory:
                    out, err, cmd = bbtools.bbduk_filter(
                        forward_in=os.path.join(args.output_folder, 'inclusion_kmers.fasta'),
                        forward_out=os.path.join(args.output_folder, 'inclusion_noplasmid.fasta'),
                        reference=args.plasmid_filtering,
                        rskip='6',
                        threads=str(args.threads),
                        returncmd=True,
                        overwrite='t', )
                    write_to_logfile(logfile=log, out=out, err=err, cmd=cmd)
                else:
                    out, err, cmd = bbtools.bbduk_filter(
                        forward_in=os.path.join(args.output_folder, 'inclusion_kmers.fasta'),
                        forward_out=os.path.join(args.output_folder, 'inclusion_noplasmid.fasta'),
                        reference=args.plasmid_filtering,
                        threads=str(args.threads),
                        returncmd=True,
                        overwrite='t')
                    write_to_logfile(logfile=log, out=out, err=err, cmd=cmd)
                # Move some sequence naming around.
                os.rename(os.path.join(args.output_folder, 'inclusion_kmers.fasta'),
                          os.path.join(args.output_folder, 'inclusion_with_plasmid.fasta'))
                os.rename(os.path.join(args.output_folder, 'inclusion_noplasmid.fasta'),
                          os.path.join(args.output_folder, 'inclusion_kmers.fasta'))
            # Now attempt to generate contiguous sequences.
            if len(glob.glob(os.path.join(args.inclusion, '*.f*a'))) > 0:
                logging.info('Generating contiguous sequences from inclusion kmers...')
                # The reference inclusion genome is the first genome returned by glob of the genomes in the inclusion
                # database folder
                ref_fasta = glob.glob(os.path.join(args.inclusion, '*.f*a'))[0]
                # Generate the coverage bedfile of the unique kmers mapped to the reference genome
                generate_bedfile(ref_fasta=ref_fasta,
                                 kmers=os.path.join(args.output_folder, 'inclusion_kmers.fasta'),
                                 output_bedfile=os.path.join(args.output_folder, 'regions_to_mask.bed'),
                                 tmpdir=os.path.join(args.output_folder, 'bedtmp'),
                                 threads=str(args.threads),
                                 logfile=log,
                                 keep=args.keep_tmpfiles)
                # Mask regions of the reference genome (with Ns) that do not have kmer coverage
                mask_fasta(input_fasta=ref_fasta,
                           output_fasta=os.path.join(args.output_folder, 'inclusion_sequence.fasta'),
                           bedfile=os.path.join(args.output_folder, 'regions_to_mask.bed'),
                           k=args.kmer_size)
                # Remove the masked regions (Ns), and create contigs from the stretches of remaining genome sequence
                remove_n(input_fasta=os.path.join(args.output_folder, 'inclusion_sequence.fasta'),
                         output_fasta=os.path.join(args.output_folder, 'sigseekr_result.fasta'),
                         k=args.kmer_size)
                # Read in any results from the output file
                with open(os.path.join(args.output_folder, 'sigseekr_result.fasta')) as f:
                    lines = f.readlines()
                # If we have file output, we have contiguous sequences long enough, so don't iterate
                if lines:
                    break
        logging.info('Could not generate signature sequences with a kmer cutoff: {cutoff}. '
                     'Trying again with cutoff: {greater}.'.format(cutoff=exclusion_cutoff,
                                                                   greater=exclusion_cutoff + 1))
        # Increment the minimum kmer cutoff value and try again
        exclusion_cutoff += 1
    # If the SigSeekr_result.fasta file has not been created, then there are no signature sequences.
    if not os.path.isfile(os.path.join(args.output_folder, 'sigseekr_result.fasta')):
        logging.info('Could not find any signature sequences in the inclusion database')
        # Don't bother trying to find primers for non-existent sequences
        args.pcr = None

    # Try to generate amplicons if desired
    if args.pcr:
        logging.info('Generating PCR info...')
        # Step 0: Create Blast DB of all exclusion genomes.
        make_all_exclusion_blast_db(exclusion_folder=args.exclusion,
                                    combined_exclusion_fasta=os.path.join(args.output_folder,
                                                                          'exclusion_combined.fasta'),
                                    logfile=log)
        for amp_size in args.amplicon_size:
            logging.info('Finding amplicons of size {}...'.format(amp_size))
            # Step 1: Go through the sigseekr_result.fasta file to find all potential amplicons based on user-specified
            # amplicon length.
            split_sequences_into_amplicons(
                input_sequence_file=os.path.join(args.output_folder, 'sigseekr_result.fasta'),
                output_amplicon_file=os.path.join(args.output_folder, 'potential_pcr_{}.fasta'.format(amp_size)),
                amplicon_length=amp_size,
                max_potential_amplicons=args.max_potential_amplicons)
            # Step 2: Blast each potential amplicon against Blast DB - keep only those that do not have any matches (for
            # now - may need to adjust this to keeping some if they have a certain e-value/length/percent id).
            ensure_amplicons_not_in_exclusion(
                exclusion_blastdb=os.path.join(args.output_folder, 'exclusion_combined.fasta'),
                potential_amplicons=os.path.join(args.output_folder, 'potential_pcr_{}.fasta'.format(amp_size)),
                confirmed_amplicons=os.path.join(args.output_folder, 'not_in_exclusion_amplicons.fasta'),
                max_potential_amplicons=args.max_potential_amplicons)
            # Step 3: Make sure potential amplicon is present in all of the inclusion genomes.
            # To do this: create blast database for each inclusion genome, and then blast each amplicon against
            # each of the inclusion genomes. Ensure that top hit is a) full length and b) pretty much identical (> 99%?)
            confirm_amplicons_in_all_inclusion_genomes(
                inclusion_fasta_dir=args.inclusion,
                potential_amplicon_file=os.path.join(args.output_folder, 'not_in_exclusion_amplicons.fasta'),
                confirmed_amplicon_file=os.path.join(args.output_folder, 'confirmed_amplicons_{}.fasta'
                                                     .format(amp_size)),
                logfile=log,
                tmpdir=os.path.join(args.output_folder, 'inclusion_pcr_tmp'),
                amplicon_size=amp_size,
                keep=args.keep_tmpfiles,
                max_potential_amplicons=args.max_potential_amplicons)
        # Now that we've generated our amplicons, we need to iterate through them and run primer3 on each, then
        # report back to the user primer pairs, amplicon sizes, and other relevant stats (melting temps, etc)
        if args.primer3:
            logging.info('Running Primer3 on potential amplicons...')
            potential_amplicon_files = glob.glob(os.path.join(args.output_folder, 'confirmed_amplicons_*.fasta'))
            run_primer3(amplicon_files=potential_amplicon_files,
                        output_file=os.path.join(args.output_folder, 'amplicons.csv'))

    if not args.keep_tmpfiles:
        logging.info('Removing unnecessary output files...')
        to_remove = glob.glob(os.path.join(args.output_folder, '*exclusion*fasta*'))
        to_remove += glob.glob(os.path.join(args.output_folder, 'unique*'))
        to_remove += glob.glob(os.path.join(args.output_folder, '*kmc*'))
        to_remove += glob.glob(os.path.join(args.output_folder, '*.bed'))
        to_remove += glob.glob(os.path.join(args.output_folder, '*sequence*'))
        for item in to_remove:
            try:
                os.remove(item)
            except FileNotFoundError:  # In case anything was already deleted, don't try to delete it twice.
                pass
    logging.info('SigSeekr run complete!')


if __name__ == '__main__':
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
    parser.add_argument('-s', '--kmer_size',
                        type=int,
                        default=31,
                        help='Kmer size used to search for sequences unique to inclusion. Default 31. No idea '
                             'how changing this affects results. TO BE INVESTIGATED.')
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
    parser.add_argument('-p', '--plasmid_filtering',
                        type=str,
                        default='NA',
                        help='To ensure unique sequences are not plasmid-borne, a FASTA-formatted database can be'
                             ' provided with this argument. Any unique kmers that are in the plasmid database will'
                             ' be filtered out.')
    parser.add_argument('-l', '--low_memory',
                        default=False,
                        action='store_true',
                        help='Activate this flag to cause plasmid filtering to use substantially less RAM (and '
                             'go faster), at the cost of some sensitivity.')
    parser.add_argument('-p3', '--primer3',
                        default=False,
                        action='store_true',
                        help='If enabled, will run primer3 on your potential amplicons and generate a list of primers '
                             'and the sizes of their products. This output will be found in a file called '
                             'amplicons.csv in the output directory specified.')
    parser.add_argument('-a', '--amplicon_size',
                        nargs='+',
                        default=[200],
                        type=int,
                        help='Desired size for PCR amplicons. Default 200. If you want to find more than one amplicon'
                             ' size, enter multiple, separated by spaces.')
    parser.add_argument('-m', '--max_potential_amplicons',
                        type=int,
                        default=200,
                        help='If inclusion sequences are very different from exclusion sequences, amplicon generation '
                             'can take forever. Set the number of potential amplicons with this option (default 200)')
    parser.add_argument('-x', '--max_memory',
                        type=int,
                        default=12,
                        help='Memory to provide (in GB) to KMC for kmerization. Default is 12')
    arguments = parser.parse_args()
    SetupLogging()
    # Check that dependencies are present, warn users if they aren't.
    dependencies = ['bbmap.sh', 'bbduk.sh', 'kmc', 'bedtools', 'samtools', 'kmc_tools', 'blastn', 'makeblastdb']
    for dependency in dependencies:
        if dependency_check(dependency) is False:
            print('WARNING: Dependency {} not found. SigSeekr may not be able to run!'.format(dependency))
    if not os.path.isdir(arguments.output_folder):
        os.makedirs(arguments.output_folder)
    main(arguments)
