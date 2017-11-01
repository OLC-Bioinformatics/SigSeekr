from argparse import ArgumentParser
import glob
from USSPpip import SigSeekr
import os
import shutil
import gzip
from All2AllMash import run_mash, read_mash
import subprocess


__author__ = 'mikeknowles'
""" This is the core script, its mission:
Retrieve genomes
This will require the user to download rMLST data
To sort rMLST results and remove closely related sequences
then to prepare the data for strain-specific probe idenification
"""


def uncompress_file(filename):
    in_gz = gzip.open(filename, 'rb')
    out = open(filename.replace('.gz', ''), 'wb')
    out.write(in_gz.read())
    out.close()


def remove_plasmid_sequences(target_file):
    """
    Also masks phage sequences, so that's cool.
    :param target_file: File you want to remove plasmid/phage sequences from.
    :return:
    """
    cmd = 'bbduk.sh ref=plasmid_database.fa in={} out={} overwrite'.format(target_file, target_file + '_1.fasta')
    print('Removing plasmid sequences...')
    print(cmd)
    subprocess.call(cmd, shell=True)
    cmd = 'bbduk.sh ref=combinedtargets.tfa in={} out={} kmask=N'.format(target_file + '_1.fasta', target_file + '_2.fasta')
    print(cmd)
    subprocess.call(cmd, shell=True)
    os.remove(target_file + '_1.fasta')
    return target_file + '_2.fasta'



def sorter(genomes, outdir, target, evalue, estop, mash_cutoff, threads):
    # if not os.path.exists(outdir + "Genomes/"):
    #     retriever(genomes, outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if not os.path.exists(outdir + "tmp/"):
        os.makedirs(outdir + "tmp/")
    # genomes = outdir + "Genomes/"
    # nontargets = glob(genomes + "*.fa*")
    to_remove = list()
    nontargets = glob.glob(genomes + '/*.f*a*')
    if len(nontargets) > 1:
        for i in range(len(nontargets)):
            if nontargets[i].endswith('.gz'):
                uncompress_file(nontargets[i])
                to_remove.append(nontargets[i].replace('.gz', ''))
                nontargets[i] = nontargets[i].replace('.gz', '')
        run_mash(genomes, threads)
        nontargets = read_mash("tmp/distances.txt", mash_cutoff)
        shutil.rmtree("tmp/")
    # If target is a folder, try to find something matching to any of the files by using a concatenated
    # sequence...not exactly what I meant to do with this. Probably need to actually go through and find a core
    # genome to get a sequence that's common to all of the files as you previously intended.
    if os.path.isdir(target):
        fasta_files = glob.glob(target + '/*.f*a*')
        for fasta in fasta_files:
            if fasta.endswith('.gz'):
                uncompress_file(fasta)
                to_remove.append(fasta.replace('.gz', ''))
                os.system('cat ' + fasta.replace('.gz', '') + ' >> ' + outdir + '/concatenated_target.fasta')
            else:
                os.system('cat ' + fasta + ' >> ' + outdir + '/concatenated_target.fasta')
            # SigSeekr(fasta, nontargets, outdir, float(evalue), float(estop), 200, 1)
        target = outdir + '/concatenated_target.fasta'
        to_remove.append(target)
        target = remove_plasmid_sequences(target)
        to_remove.append(target)
    else:
        if '.gz' in target:
            uncompress_file(target)
            to_remove.append(target.replace('.gz', ''))
            target = target.replace('.gz', '')
            target = remove_plasmid_sequences(target)
            to_remove.append(target)
        else:
            target = remove_plasmid_sequences(target)
            to_remove.append(target)
    # else:
    SigSeekr(target, nontargets, outdir, float(evalue), float(estop), 200, 1)
    for item in to_remove:
        os.remove(item)

if __name__ == '__main__':
    import multiprocessing
    num_cpus = multiprocessing.cpu_count()
    # Parser for arguments
    parser = ArgumentParser(description='Find Universal Strain-Specifc Probes')
    parser.add_argument('--version', action='version', version='%(prog)s v0.5')
    parser.add_argument('-o', '--output', required=True, help='Specify output directory')
    parser.add_argument('-i', '--input', required=True, help='Specify input genome fasta folder')
    parser.add_argument('-t', '--target', required=True, help='Specify target genome or folder')
    parser.add_argument('-e', '--evalue', type=float, default=1e-40, help='Specify elimination E-value lower limit (default 1e-50)')
    parser.add_argument('-s', '--estop', type=float, default=1e-90, help='Specify the upper E-value limit (default 1e-90)')
    parser.add_argument('-c', '--mash_cutoff', type=float, default=0.0002, help='Cutoff value to use for genome'
                                                                                ' elimination. Must be a float between'
                                                                                '0 and 1. Default is 0.0002, higher values'
                                                                                'eliminate more genomes.')
    parser.add_argument('-n', '--num_threads', type=int, default=num_cpus, help='Number of threads to run analysis on.'
                                                                                'Defaults to number of CPUs on system.')
    # parser.add_argument('-t', '--target', required=True, help='Specify target genome or folder')
    args = vars(parser.parse_args())

    sorter(
           os.path.join(os.path.abspath(args['input']), ""),
           os.path.join(os.path.abspath(args['output']), ""),
           os.path.abspath(args['target']), args['evalue'],
       args['estop'], args['mash_cutoff'], args['num_threads'])