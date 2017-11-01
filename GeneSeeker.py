__author__ = 'mikeknowles'
""" Includes threading found in examples:
http://www.troyfawkes.com/learn-python-multithreading-queues-basics/
http://www.ibm.com/developerworks/aix/library/au-threadingpython/
https://docs.python.org/2/library/threading.html
Revised with speed improvements
"""
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from threading import Thread
from queue import Queue
from collections import defaultdict
from io import StringIO
import subprocess, os, glob, time, sys, shlex, re, threading, json, mmap
count = 0


def dotter():
    global count
    # if count <= 80:
    sys.stdout.write('.')
        # count += 1
    # else:
    #     sys.stdout.write('\n[%s].' % (time.strftime("%H:%M:%S")))
    #     count = 0


def makeblastdb(dqueue):
    while True:  # while daemon
        fastapath = dqueue.get() # grabs fastapath from dqueue
        #db = fastapath  # remove the file extension for easier future globing
        db = fastapath.replace(".fasta","")
        nhr = "%s.nhr" % db  # add nhr for searching
        FNULL = open(os.devnull, 'w')  # define /dev/null
        if not os.path.isfile(str(nhr)):  # if check for already existing dbs
            subprocess.Popen(shlex.split("makeblastdb -in %s -dbtype nucl -out %s" % (fastapath, db)), stderr=FNULL, stdout=FNULL)
            # make blastdb
            dotter()
        dqueue.task_done() # signals to dqueue job is done

# Declare queues, list, and dict
dqueue = Queue()
blastqueue = Queue()
parsequeue = Queue()
testqueue = Queue()
plusqueue = Queue()
plusdict = {}
genedict = defaultdict(list)
blastpath = []
# threadlock, useful for overcoming GIL
threadlock = threading.Lock()

def makedbthreads(fastas):
    ''' Setup and create threads for class'''
    for i in range(12):
        threads = Thread(target=makeblastdb, args=(dqueue,))
        threads.setDaemon(True)
        threads.start()
    for fasta in fastas:
        dqueue.put(fasta)
    dqueue.join()  # wait on the dqueue until everything has been processed

def xmlout (fasta, genome):
    path = re.search('(.+)\/(.+)\/(.+?)$', genome)
    gene = fasta.split('/')[-1]  # split file from path, could use os.path.split
    genename = gene.split('.')[0]
    genomename = path.group(3).split('.')[0]
    out = "%s/tmp/%s.%s.xml" % (path.group(1), genomename, genename)  # Changed from dictionary to tuple
    return path, gene, genename, genomename, out

def blastparse(blast_handle, genome, gene):
    global plusdict
    records = NCBIXML.parse(blast_handle)   # Open record from memory-mapped file
    dotter()
    for record in records:  # This process is just to retrieve HSPs from xml files
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                threadlock.acquire()  # precaution
                # if hsp.identities == alignment.length:  # if the length of the match matches the legth of the sequence
                #     # if genome not in plusdict:  # add genomes in plusdict
                #     #     plusdict[genome] = defaultdict(list)
                #     # if gene not in plusdict[genome]:  # add genes to plus dict
                #     #     plusdict[genome][gene] = []
                if plusdict[genome][gene] == [] and hsp.identities == alignment.length:
                    # If there is only one good match then apply allele number
                    plusdict[genome][gene].append(alignment.title.split('_')[-1])
                elif hsp.identities == alignment.length:
                    # If there is multiple matches then added them in a string
                    plusdict[genome][gene].append(alignment.title.split('_')[-1])
                    plusdict[genome][gene].sort()
                else:
                    # or add the
                    plusdict[genome][gene].append('%s (%s/%s)' % (alignment.title.split('_')[-1],
                                                                  hsp.identities,
                                                                  alignment.length))

                threadlock.release()  # precaution for populate dictionary with GIL



class runblast(threading.Thread):
    def __init__(self, blastqueue):
        self.blastqueue = blastqueue
        threading.Thread.__init__(self)
    def run(self):
        while True:
            global blastpath, plusdict  # global varibles, might be a better way to pipe information
            genome, fasta, blastexist = self.blastqueue.get()  # retrieve variables from queue
            path, gene, genename, genomename, out = xmlout(fasta, genome)  # retrieve from string splitter
            threadlock.acquire()
            blastpath.append((out, path.group(3), gene, genename,))  # tuple-list
            try:
                plusdict[genome][genename] = []
            except KeyError:
                plusdict[genome] = defaultdict(list)
                plusdict[genome][genename] = []
            threadlock.release()
            if not os.path.isfile(out):
                dotter()
                file = open(out, 'w')
                # for perc in range(100, 99, -1):  # try to get allele type at varying ident
                db = fasta.split('.')[0]
                blastn = NcbiblastnCommandline(query=genome, db=db, evalue=1e-40, outfmt=5, perc_identity=100)
                stdout, stderr = blastn()
                if stdout.find('Hsp') != -1:
                    blast_handle = StringIO(stdout)  # Convert string to IO object for use in SearchIO using StringIO
                    blastparse(blast_handle, genome, genename)  # parse the data already in memory
                    file.write(stdout)  # write the result
                        # break
                file.close()
            elif os.path.getsize(out) != 0:  # if blast files already exist
                handle = open(out)
                mm = mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ)
                parsequeue.put((out, genome, genename, mm))
            parsequeue.join()
            if not any(blastpath):
                print(out)
            self.blastqueue.task_done()


def blastnthreads(fastas, genomes):
    '''Setup and create  threads for blastn and xml path'''
    blastexist = {}

    for genome in genomes:
        for fasta in fastas:
            blastqueue.put((genome, fasta, blastexist))
        blastqueue.join()


class multiparser(threading.Thread): # Had to convert this to a class to integrate threading lock
    def __init__(self, parsequeue):
        self.parsequeue = parsequeue
        threading.Thread.__init__(self)

    def run(self):
        while True:  # General Loop
            global plusdict, genedict  # Import global elements to populate, there may be a better way to do this
            xml, genome, gene, mm = self.parsequeue.get()  # Retrieve dara from queue
            blastparse(mm, genome, gene)
            mm.close()
            # TODO: Add genefinder-like functionality here using a queue
            self.parsequeue.task_done()

def parsethreader(blastpath):
    global plusdict
    for pathstr in blastpath:
        xml = pathstr[0]
        if os.path.getsize(xml) != 0:
            handle = open(xml, 'r')
            mm = mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ)
            # time.sleep(0.05) # Previously used to combat open file error
            handle.close()
            parsequeue.put((xml, pathstr[1], pathstr[3], mm))
            parsequeue.join()


def blaster(markers, strains, out, name):
    '''
    The blaster function is the stack manager of the module
    markers are the the target fasta folder that with be db'd and BLAST'd against strains folder
    out is the working directory where the blastxml folder will be placed
    name is the partial title of the csv output
    ALL PATHS REQUIRE TRAILING SLASHES!!!
    '''
    # TODO: Keep file extension when making dictionary and splice end of
    start = time.time()
    global count, genedict, blastpath, plusdict
    plusdict = {}
    genedict = defaultdict(list)
    blastpath = []

    #retrieve markers from input
    genes = glob.glob(markers + "*.fas")
    for i in range(len(strains)):
        threads = runblast(blastqueue)
        threads.setDaemon(True)
        parthreads = multiparser(parsequeue)
        parthreads.setDaemon(True)
        threads.start()
        parthreads.start()
    #retrieve genomes from input
    if os.path.isdir(strains):
        genomes = glob.glob(strains + "*.f*")
        print('[%s] GeneSeekr input is path with %s genomes' % (time.strftime("%H:%M:%S"), len(genomes)))
    elif os.path.isfile(strains):
        genomes = [strains,]
        strains = os.path.split(strains)[0]
        print('GeneSeeker input is a single file \n%s' % genomes)
    else:
        print("The variable \"--genomes\" is not a folder or file")
        return
    sys.stdout.write("[%s] Creating necessary databases for BLAST" % (time.strftime("%H:%M:%S")))
    #push markers to threads
    makedbthreads(genes)
    print("\n[%s] BLAST database(s) created" % (time.strftime("%H:%M:%S")))
    if os.path.isfile('%s/%s_blastxmldict.json' % (strains, name)):
        print("[%s] Loading BLAST data from file" % (time.strftime("%H:%M:%S")))
        blastpath = json.load(open('%s/%s_blastxmldict.json' % (strains, name)))
        print("[%s] Now parsing BLAST database searches" % (time.strftime("%H:%M:%S")))
        sys.stdout.write('[%s] ' % (time.strftime("%H:%M:%S")))
        parsethreader(blastpath)
    else:
        print("[%s] Now performing and parsing BLAST database searches" % (time.strftime("%H:%M:%S")))
        sys.stdout.write('[%s]' % (time.strftime("%H:%M:%S")))
        # make blastn threads and retrieve xml file locations
        blastnthreads(genes, genomes)
        print('\n')
        json.dump(blastpath, open('%s/%s_blastxmldict.json' % (strains, name), 'w'), sort_keys=True, indent=4, separators=(',', ': '))
    csvheader = 'Strain'
    row = ""
    rowcount = 0
    # create csv file with strings
    for genomerow in plusdict:
        row += "\n" + genomerow.split('/')[-1].split('.')[0]
        rowcount += 1

        for generow in sorted(genes):
            genename = generow[-14:-4]
            if rowcount <= 1:
                csvheader += ', ' + genename
            if genename in plusdict[genomerow]:
                alleleNum = ""
                for allele in plusdict[genomerow][genename]:
                    alleleNum += str(allele) + ' '
                row += ',' + alleleNum
            else:
                row += ',N'
    with open("%s%s_results_%s.csv" % (out, name, time.strftime("%Y.%m.%d.%H.%M.%S")), 'wb') as csvfile:
        csvfile.write(csvheader)
        csvfile.write(row)
    end = time.time() - start
    print("\n[%s] Elapsed time for rMLST is %ss with %ss per genome" % (time.strftime("%H:%M:%S"), end, end/float(len(genomes))))
    return plusdict
