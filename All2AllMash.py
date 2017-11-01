# Takes some stuff in, mashes everything against everything.

import os


def run_mash(folder, num_threads):
    # Find all the fasta files in the user-specified folder.
    # Sketch the files, using number of threads specified. Put the output into tmp.
    if not os.path.isdir(os.path.join(os.getcwd(),"tmp")):
        os.makedirs("tmp")
    cmd = "mash sketch -p " + str(num_threads) + " -o tmp/reference.msh " + folder + "*.f*a"
    os.system(cmd)
    # Now that things are sketched, do the all to all. Results also go into tmp.
    cmd = "mash dist -p " + str(num_threads) + " tmp/reference.msh " + folder + "*.f*a" + " > tmp/distances.txt"
    os.system(cmd)


def read_mash(distance_file, cutoff):
    """
    :param distance_file: Distances file created by mash.
    :param cutoff: Mash distance that marks genomes that are close together.
    :return: to_keep: list of genomes that are far enough apart to be worthwhile.
    """
    to_remove = list()
    to_keep = list()
    f = open(distance_file)
    lines = f.readlines()
    f.close()

    # (I.e. if A and B are close, and B and C are close, but A and C are far, how do we eliminate only B consistently).
    # This seems to be working now, so that's good. We'll get this integrated into SigSeekr in the future.
    for line in lines:
        x = line.split()
        if x[0] != x[1]:
            if float(x[2]) <= cutoff:
                if x[1] not in to_remove and x[0] not in to_remove:
                    to_remove.append(x[1])

    for line in lines:
        x = line.split()
        if x[0] not in to_remove and x[0] not in to_keep:
            to_keep.append(x[0])

    return to_keep

"""
if __name__ == "__main__":
    import argparse
    import time
    import multiprocessing
    start = time.time()
    cpu_count = multiprocessing.cpu_count()
    parser = argparse.ArgumentParser()
    parser.add_argument("fastafolder", help="Folder that contains fasta files you want to mash. " 
                        "Extension must be .fasta")
    arguments = parser.parse_args()
    run_mash(arguments.fastafolder, cpu_count)
    read_mash("distances.txt", 0.0002)
    end = time.time()
    m, s = divmod(end - start, 60)
    h, m = divmod(m, 60)
    print("Finished mashing in %d:%02d:%02d " % (h, m, s))
"""
