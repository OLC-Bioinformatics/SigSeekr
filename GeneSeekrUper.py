#! /usr/env/python
__author__ = 'mikeknowles, akoziol'
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
from glob import glob
import subprocess, os, time, sys, shlex, re, threading, json, mmap, errno
from argparse import ArgumentParser

# parser = ArgumentParser(description='Performs blast analyses to determine presence of custom targets')
# parser.add_argument('-p', '--path', required=False,
#                     default='/home/blais/PycharmProjects/pythonGeneSeekr/',
#                     help='Specify path for custom folder locations')
# parser.add_argument('-c', '--cutoff', required=False, default=0.8,
#                     help='The identity cutoff value for BLAST matches. Default is 0.8')
# parser.add_argument('-s', '--sequencePath', required=False,
#                     default='/home/blais/PycharmProjects/pythonGeneSeekr/sequences',
#                     help='The location of the query sequence files')
# parser.add_argument('-t', '--targetPath', required=False,
#                     default='/home/blais/PycharmProjects/pythonGeneSeekr/Organism',
#                     help='The location of the target files')
#
# # Get the arguments into a list
# args = vars(parser.parse_args())

# Define variables from the arguments - there may be a more streamlined way to do this
# Add trailing slashes to the path variables to ensure consistent formatting (os.path.join)
# path = os.path.join(args['path'], "")
# cutoff = float(args['cutoff'])
# sequencePath = os.path.join(args['sequencePath'], "")
# targetPath = os.path.join(args['targetPath'], "")

def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
        # os.chmod(inPath, 0777)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def make_dict():
    """Makes Perl-style dictionaries"""
    return defaultdict(make_dict)

# Initialise the count used in the dotter function
count = 0


def dotter():
    """Prints formatted time to stdout at the start of a line, as well as a "."
    whenever the length of the line is equal or lesser than 80 "." long"""
    # Use a global variable
    global count
    if count <= 80:
        sys.stdout.write('.')
        count += 1
    else:
        sys.stdout.write('\n[%s] .' % (time.strftime("%H:%M:%S")))
        count = 1


def makeblastdb(dqueue):
    while True:  # while daemon
        fastapath = dqueue.get()  # grabs fastapath from dqueue
        # remove the path and the file extension for easier future globbing
        db = fastapath.split(".")[0]
        nhr = "%s.nhr" % db  # add nhr for searching
        # print nhr
        FNULL = open(os.devnull, 'w')  # define /dev/null
        if not os.path.isfile(str(nhr)):  # if check for already existing dbs
            subprocess.Popen(shlex.split("makeblastdb -in %s -dbtype nucl -out %s" % (fastapath, db)))
            # make blastdb
            dotter()
        dqueue.task_done()  # signals to dqueue job is done
        sys.exit()  # necessary

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
    """Setup and create threads for class"""
    # Create and start threads for each fasta file in the list
    for i in range(len(fastas)):
        # Send the threads to makeblastdb
        threads = Thread(target=makeblastdb, args=(dqueue,))
        # Set the daemon to true - something to do with thread management
        threads.setDaemon(True)
        # Start the threading
        threads.start()
    for fasta in fastas:
        # Add the fasta file to the queue
        dqueue.put(fasta)
    dqueue.join()  # wait on the dqueue until everything has been processed


def xmlout(fasta, genome):
    """Parses variables from supplied tuples? dictionaries?"""
    path = genome.split("/")
    # path = re.search('(.+)\/(.+)\/(.+?)$', genome)
    # print path
    # print path.group(2)
    gene = fasta.split('/')[-1]  # split file from path, could use os.path.split
    genename = gene.split('.')[0]
    genomename = path[-1].split('.')[0]
    # print genomename
    # Create the out variable containing the path and name of BLAST output file
    tmpPath = path[0]
    # tmpPath = "%s/%s" % (path.group(1), path.group(2))
    make_path("%s/tmp" % tmpPath)
    out = "%s/tmp/%s.%s.xml" % (tmpPath, genomename, genename)  # Changed from dictionary to tuple
    # Return the parsed variables
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
                genename = alignment.title.split(" ")[1].split('_')
                gene = genename[0]
                if gene not in plusdict[genome]:
                    plusdict[genome][gene] = []
                if plusdict[genome][gene] == [] and hsp.identities == alignment.length:
                    # If there is only one good match then apply allele number
                    plusdict[genome][gene].append(genename[-1])
                elif hsp.identities == alignment.length:
                    # If there is multiple matches then added them in a string
                    plusdict[genome][gene].append(genename[-1])
                    plusdict[genome][gene].sort()
                else:
                    # or add the
                    plusdict[genome][gene].append('%s (%s/%s)' % (genename[-1],
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
            genome, fasta, blastexist, analysisType, cutoff = self.blastqueue.get()  # retrieve variables from queue
            path, gene, genename, genomename, out = xmlout(fasta, genome)  # retrieve from string splitter
            #Precaution
            threadlock.acquire()
            # Add the appropriate variables to blast path
            blastpath.append((out, path[-1], gene, genename,))  # tuple-list
            if genome not in plusdict:
                plusdict[genome] = {}
            # try:
            #     plusdict[genome][genename] = []
            # except KeyError:
            #     plusdict[genome] = {}
            #     plusdict[genome][genename] = []
            threadlock.release()
            # Checks to see if this BLAST search has previously been performed
            if not os.path.isfile(out):
                # Print a dot for each gene, genome combination processed
                dotter()
                # Open the output file for writing
                file = open(out, 'w')
                # Run the BioPython BLASTn module with the genome as query, fasta(target gene) as db,
                # a mild evalue of 0.1, and XML formatted output
                # Removed perc_identity=percentIdentity from the call, as this allows more flexibility for parsing files
                # with different cutoff values - if I want to re-analyse a search with a lower cutoff, I won't have to
                # re-perform the BLAST search each time
                db = fasta.split(".")[0]
                blastn = NcbiblastnCommandline(query=genome, db=db, evalue=0.1, outfmt=5, out=genome + ".csv")
                # Note that there is no output file specified -  the search results are currently stored in stdout
                stdout, stderr = blastn()
                # Search stdout for matches - if the term Hsp appears (the .find function will NOT
                # return -1), a match has been found, and stdout is written to file
                if stdout.find('Hsp') != -1:
                    blast_handle = StringIO(stdout)  # Convert string to IO object for use in SearchIO using StringIO
                    blastparse(blast_handle, genome, genename)  # parse the data already in memory
                    file.write(stdout)  # write the result
                # Close the file
                file.close()
            # If the BLAST results file already exists is not empty, then parse the results
            elif os.path.getsize(out) != 0:
                # Open the file
                handle = open(out)
                # Read the file into memory
                mm = mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ)
                # Parse the file in a multithreaded manner
                parsequeue.put((out, genome, genename, mm, analysisType, cutoff))
            # Join all the threads
            parsequeue.join()
            # Error catching?
            if not any(blastpath):
                print(out)
            self.blastqueue.task_done()


def blastnthreads(fastas, genomes, analysisType, cutoff):
    """Setup and create  threads for blastn and xml path"""
    blastexist = {}
    # Create threads for each gene, genome combination
    for fasta in fastas:
        for genome in genomes:

            # Add the appropriate variables to the threads
            blastqueue.put((genome, fasta, blastexist, analysisType, cutoff))
        blastqueue.join()


class multiparser(threading.Thread): # Had to convert this to a class to integrate threading lock
    def __init__(self, parsequeue):
        self.parsequeue = parsequeue
        threading.Thread.__init__(self)

    def run(self):
        while True:  # General Loop
            global plusdict, genedict  # Import global elements to populate, there may be a better way to do this
            xml, genome, gene, mm, analysisType, cutoff = self.parsequeue.get()  # Retrieve dara from queue
            blastparse(mm, genome, gene)
            mm.close()
            self.parsequeue.task_done()


def organismChooser(path,targetPath, name):
    """Allows the user to choose which organism to be used in the analyses"""
    # Initialise a count variable to be used in extracting the desired entry from a list of organisms
    count = 0
    # Check to see if the supplied targetPath has .fa files - if it does, then the default directory structure is probably
    # not being followed, so the target files will be in targetPath
    foldertest = glob("%s/*.fa*" % targetPath)
    if foldertest:
        # Set the required variables as necessary
        # queryGenes are presumably the genes found by foldertest
        queryGenes = foldertest
        # There are likely not going to be qualityGenes included in a custom analysis
        qualityGenes = []
        # Organism name is simply the name of the folder containing the targets
        organismName = name
    else:
        # Get a list of the organisms in the (default) Organism subfolder
        orgList = glob("%sOrganism/*" % path)
        # Iterate through the sorted list
        for folder in sorted(orgList):
            # Ensure that folder is, in actuality, a folder
            if os.path.isdir(folder):
                # Print out the folder names and the count
                print("[%s]: %s" % (count, os.path.split(folder)[1]))
                count += 1
        # Get the user input - the number entered corresponds to the list index
        response = input("Please select an organism: ")
        # Get the organism path into a variable
        organism = sorted(orgList)[int(response)]
        organismName = os.path.split(organism)[1]
        # Put the query and quality genes into lists
        queryGenes = glob("%s/query_genes/*.fa" % organism)
        qualityGenes = glob("%s/qualityTest/*.tfa" % organism)
    return queryGenes, qualityGenes, organismName


def blaster(path, cutoff, sequencePath, targetPath, name):
    """
    The blaster function is the stack manager of the module
    markers are the the target fasta folder that with be db'd and BLAST'd against strains folder
    out is the working directory where the blastxml folder will be placed
    name is the partial title of the csv output
    ALL PATHS REQUIRE TRAILING SLASHES!!!
    """
    # Time is used to calculate length of the analyses
    start = time.time()
    # Import global variables
    global count, genedict, blastpath, plusdict
    # Initialise genedict
    genedict = defaultdict(list)
    blastpath = []
    # Run organism chooser to allow the user to choose which databases to use
    # returns the organism name, and lists of
    # queryGenes, qualityGenes, organismName = organismChooser(path, targetPath, name)
    # Set the required variables as necessary
    if os.path.isfile(targetPath):
        foldertest = [targetPath,]
    else:
        foldertest = glob(targetPath + "/*.fa*")
    # queryGenes are presumably the genes found by foldertest
    queryGenes = foldertest
    # There are likely not going to be qualityGenes included in a custom analysis
    qualityGenes = []
    # Organism name is simply the name of the folder containing the targets
    organismName = name
    # Get the genome files into a list - note that they must be in the "sequences" subfolder of the path,
    # and the must have a file extension beginning with ".fa"
    if os.path.isdir(sequencePath):
        strains = glob("%s*.fa*" % sequencePath)
        print('[%s] GeneSeekr input is path with %s genomes' % (time.strftime("%H:%M:%S"), len(strains)))
    elif os.path.isfile(sequencePath):
        strains = [sequencePath,]
        singlestrain = os.path.split(sequencePath)[1]
        print('GeneSeeker input is a single file \n%s' % singlestrain)
    else:
        print("The variable \"--genomes\" is not a folder or file")
        return
    # Create the threads for the BLAST analysis
    for i in range(len(strains)):
        threads = runblast(blastqueue)
        threads.setDaemon(True)
        parthreads = multiparser(parsequeue)
        parthreads.setDaemon(True)
        threads.start()
        parthreads.start()
    sys.stdout.write("[%s] Creating necessary databases for BLAST" % (time.strftime("%H:%M:%S")))
    # Push targets to threads
    makedbthreads(queryGenes)
    # Quality test genes are optional, so only run the qualityGenes if the folder exists
    if qualityGenes:
        makedbthreads(qualityGenes)
    print("\n[%s] BLAST database(s) created" % (time.strftime("%H:%M:%S")))
    print("[%s] Now performing and parsing BLAST database searches" % (time.strftime("%H:%M:%S")))
    sys.stdout.write('[%s] ' % (time.strftime("%H:%M:%S")))
    # Make blastn threads and retrieve xml file locations
    blastnthreads(queryGenes, strains, "query", cutoff)
    # qualityGenes optional
    if qualityGenes:
        blastnthreads(qualityGenes, strains, "quality", cutoff)
    # Initialise types dictionary
    types = {}
    # Populate types
    types["query"] = queryGenes
    if qualityGenes:
        types["quality"] = qualityGenes
    # csvheader = ''
    # Loop through the analysis types, and make outputs as required
    # for analysisType in types:
    #     # Initialise variables
    #     row = ""
    #     rowcount = 0
    #     # plusdict contains all results - format: plusdict[genomeName][gene][analysisType] = totalPercentID
    #     for genomerow in plusdict:
    #         # The first cell contains the word "Strain"
    #         csvheader = 'Strain'
    #         # Append the genome names taken from genomerow
    #         row += "\n" + genomerow.split('/')[-1].split('.')[0]
    #         # Increment rowcount
    #         rowcount += 1
    #         # For each gene in the appropriate analysis type
    #         for generow in sorted(types[analysisType]):
    #             # Extract the gene name from generow
    #             genename = os.path.basename(generow).split('.')[0]
    #             # If the gene name is not already in the csvheader variable
    #             if genename not in csvheader:
    #                 # Append the gene name in a comma-separated format
    #                 csvheader += ',' + genename
    #                 # Format the results - if the percent ID is less than 100
    #                 identity = plusdict[genomerow][genename][analysisType]
    #                 if cutoff * 100 < float(identity) < 100:
    #                 #  Append the percent ID (and a "%") to the row variable
    #                     row += ',' + str(plusdict[genomerow][genename][analysisType]) + "%"
    #                 # Otherwise, add a "+" to represent a 100% match
    #                 elif identity == 0:
    #                     row += ',N'
    #                 # If the analysisType is 0 in the dictionary, then there were no matches.
    #                 # This is shown by an 'N'
    #                 else:
    #                     row += ',+'
    #     # Open the csv report in the appropriate location - add the organism name and the date to keep reports unique
    #     make_path("%sreports" % path)
    #     with open("%sreports/%s_%s_results_%s.csv" % (path, organismName, analysisType, time.strftime("%Y.%m.%d.%H.%M.%S")), 'wb') as csvfile:
    #         # Write the header and the rows
    #         csvfile.write(csvheader)
    #         csvfile.write(row)
    # # Calculate the elapsed time
    end = time.time() - start
    # # Friendly exit statement
    print("\n[%s] Elapsed time for GeneSeeking is %.2f seconds with %.2f seconds per genome" \
          % (time.strftime("%H:%M:%S"), end, end/float(len(strains))))
    return plusdict

# Run the blaster function
# blaster(path, cutoff, sequencePath, targetPath)
