from biotools import kmc
import sys
import os
import shutil

parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0, parentdir)
from sigseekr.sigseekr import *


def test_replace_by_index():
    seq = 'ATGTTACGTAACCTGTACCATGTACACAGGT'
    seq = replace_by_index('1:5', seq)
    assert seq == 'ANNNNACGTAACCTGTACCATGTACACAGGT'


def test_n_removal_short():
    remove_n('tests/n_removal.fasta', 'tests/n_result.fasta')
    with open('tests/n_result.fasta') as f:
        lines = f.readlines()
    assert lines == []
    os.remove('tests/n_result.fasta')


def test_n_removal_long():
    remove_n('tests/n_removal_long.fasta', 'tests/n_result.fasta')
    with open('tests/n_result.fasta') as f:
        lines = f.readlines()
    assert lines == ['>contig1_sequence1\n', 'TGCATCGAATATGCATGCAATGCATGCAATGCATGCAGCATGCAAAAAAAAAA\n']


def test_1_fastqs():
    assert find_paired_reads('tests/fake_fastqs/', forward_id='_1',
                             reverse_id='_2') == [['tests/fake_fastqs/test_1.fastq.gz',
                                                   'tests/fake_fastqs/test_2.fastq.gz']]


def test_unpaired_fastqs():
    unpaired = find_unpaired_reads('tests/fake_fastqs')
    assert 'tests/fake_fastqs/test_R1.fastq.gz' in unpaired
    assert 'tests/fake_fastqs/test_alone.fastq.gz' in unpaired


def test_inclusion_fastas():
    if not os.path.isdir('tests/tmp'):
        os.makedirs('tests/tmp')
    make_inclusion_kmerdb('tests/fasta_only', 'tests/tmp/testdb', maxmem='1')
    kmc.dump('tests/tmp/testdb', 'tests/tmp/kmers')
    with open('tests/tmp/kmers') as f:
        lines = f.readlines()
    assert lines == ['ATGTAGCATGTACGTACGTAGTCATGGGATA	1\n']
    shutil.rmtree('tests/tmp')


def test_exclusion_fastas():
    if not os.path.isdir('tests/tmp'):
        os.makedirs('tests/tmp')
    make_exclusion_kmerdb('tests/fasta_only', 'tests/tmp/testdb', maxmem='1')
    kmc.dump('tests/tmp/testdb', 'tests/tmp/kmers')
    with open('tests/tmp/kmers') as f:
        lines = f.readlines()
    assert lines == ['ATGTAGCATGTACGTACGTAGTCATGGGATA	2\n',
                     'CCCTATCCCATGACTACGTACGTACATGCTA	1\n',
                     'CCTATCCCATGACTACGTACGTACATGCTAC	1\n',
                     'CTATCCCATGACTACGTACGTACATGCTACA	1\n',
                     'GGGTATCCCATGACTACGTACGTACATGCTA	1\n',
                     'GGTATCCCATGACTACGTACGTACATGCTAC	1\n',
                     'GTATCCCATGACTACGTACGTACATGCTACA	1\n']
    shutil.rmtree('tests/tmp')


def test_inclusion_fastqs():
    if not os.path.isdir('tests/tmp'):
        os.makedirs('tests/tmp')
    make_inclusion_kmerdb('tests/fastq_only', 'tests/tmp/testdb', maxmem='1')
    kmc.dump('tests/tmp/testdb', 'tests/tmp/kmers')
    with open('tests/tmp/kmers') as f:
        lines = f.readlines()
    assert lines == ['ATGTAGCATGTACGTACGTAGTCATGGGATA	2\n']
    shutil.rmtree('tests/tmp')


def test_exclusion_fastqs():
    if not os.path.isdir('tests/tmp'):
        os.makedirs('tests/tmp')
    make_exclusion_kmerdb('tests/fastq_only', 'tests/tmp/testdb', maxmem='1')
    kmc.dump('tests/tmp/testdb', 'tests/tmp/kmers')
    with open('tests/tmp/kmers') as f:
        lines = f.readlines()
    assert lines == ['AAATATCCCATGACTACGTACGTACATGCTA	2\n',
                     'AATATCCCATGACTACGTACGTACATGCTAC	2\n',
                     'ATATCCCATGACTACGTACGTACATGCTACA	2\n',
                     'ATGTAGCATGTACGTACGTAGTCATGGGATA	4\n',
                     'GTAGCATGTACGTACGTAGTCATGGGATAAA	2\n',
                     'TAGCATGTACGTACGTAGTCATGGGATAAAA	2\n',
                     'TGTAGCATGTACGTACGTAGTCATGGGATAA	2\n']
    shutil.rmtree('tests/tmp')


def test_inclusion_both():
    if not os.path.isdir('tests/tmp'):
        os.makedirs('tests/tmp')
    make_inclusion_kmerdb('tests/fasta_and_fastq', 'tests/tmp/testdb', maxmem='1')
    kmc.dump('tests/tmp/testdb', 'tests/tmp/kmers')
    with open('tests/tmp/kmers') as f:
        lines = f.readlines()
    assert lines == ['ATGTAGCATGTACGTACGTAGTCATGGGATA	1\n']
    shutil.rmtree('tests/tmp')


def test_exclusion_both():
    if not os.path.isdir('tests/tmp'):
        os.makedirs('tests/tmp')
    make_exclusion_kmerdb('tests/fasta_and_fastq', 'tests/tmp/testdb', maxmem='1')
    kmc.dump('tests/tmp/testdb', 'tests/tmp/kmers')
    with open('tests/tmp/kmers') as f:
        lines = f.readlines()
    assert lines == ['AAATATCCCATGACTACGTACGTACATGCTA	2\n',
                     'AATATCCCATGACTACGTACGTACATGCTAC	2\n',
                     'ATATCCCATGACTACGTACGTACATGCTACA	2\n',
                     'ATGTAGCATGTACGTACGTAGTCATGGGATA	6\n',
                     'CCCTATCCCATGACTACGTACGTACATGCTA	1\n',
                     'CCTATCCCATGACTACGTACGTACATGCTAC	1\n',
                     'CTATCCCATGACTACGTACGTACATGCTACA	1\n',
                     'GGGTATCCCATGACTACGTACGTACATGCTA	1\n',
                     'GGTATCCCATGACTACGTACGTACATGCTAC	1\n',
                     'GTAGCATGTACGTACGTAGTCATGGGATAAA	2\n',
                     'GTATCCCATGACTACGTACGTACATGCTACA	1\n',
                     'TAGCATGTACGTACGTAGTCATGGGATAAAA	2\n',
                     'TGTAGCATGTACGTACGTAGTCATGGGATAA	2\n']
    shutil.rmtree('tests/tmp')


def test_kmers_to_fasta():
    kmers_to_fasta('tests/kmers.txt', 'tests/output.fasta')
    with open('tests/output.fasta') as f:
        lines = f.readlines()
    assert lines == ['>kmer1\n',
                     'ATGTAGCATGTACGTACGTAGTCATGGGATA\n']
    os.remove('tests/output.fasta')

# The following test works locally, but for some reason Travis hates it. May be related to versioning of software.
# Will have to investigate to see if it's actually that causing breaks (bedtools/bbmap/samtools are the ones involved)
"""
def test_bedfile_generation():
    generate_bedfile('tests/bed_creation/reference.fasta', 'tests/bed_creation/query.fasta', 'tests/bed_creation/test.bed')
    with open('tests/bed_creation/test.bed') as f:
        lines = f.readlines()
    assert lines == ['sequence	0	61	0\n',
                     'sequence	61	101	1\n',
                     'sequence	101	165	0\n']
    os.remove('tests/bed_creation/test.bed')
"""

def test_fasta_mask():
    mask_fasta('tests/fasta_mask/reference.fasta', 'tests/fasta_mask/out.fasta', 'tests/fasta_mask/test.bed')
    with open('tests/fasta_mask/out.fasta') as f:
        lines = f.readlines()
    assert lines == ['>sequence\n',
                     'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTAGCTACGATCGATCATCATCGAT'
                     'CATCGATCACGATCATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n']
    os.remove('tests/fasta_mask/out.fasta')


def test_amplicon():
    find_primer_distances('tests/amplicon_test/kmers.fasta', 'tests/amplicon_test/reference.fasta',
                          'tests/amplicon_test', min_amplicon_size=100,
                          inclusion_dir='tests/amplicon_test/inclusion')  # TODO: Add inclusion dir so that this works again
    with open('tests/amplicon_test/amplicons.csv') as f:
        lines = f.readlines()
    # assert lines == ['Sequence1,Sequence2,Amplicon_Size\n',
    #                  'AGCTAGTCACTACAGCTACGATAGCTAGCAT,GGCACCATCGACACGACACCACACACGACAG,123\n',
    #                  'GGCACCATCGACACGACACCACACACGACAG,AGCTAGTCACTACAGCTACGATAGCTAGCAT,123\n']
    assert 'AGCTAGTCACTACAGCTACGATAGCTAGCAT,GGCACCATCGACACGACACCACACACGACAG,123\n' in lines
    assert 'GGCACCATCGACACGACACCACACACGACAG,AGCTAGTCACTACAGCTACGATAGCTAGCAT,123\n' in lines
    assert len(lines) == 3
    os.remove('tests/amplicon_test/amplicons.csv')
