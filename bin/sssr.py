#!/usr/bin/python3

import os
import itertools
import multiprocessing as mp
import shlex
import subprocess
import sys

EXE = './release/src/main'


def init(l):
    global lock
    lock = l


def main(files):
    # Filter out any files with 'cycles' as a substring in the base filename
    files = [f for f in files if 'cycles' not in os.path.basename(f)]

    notCSVFiles      = []
    dneFiles         = []
    notGraphFiles    = []
    noOverwriteFiles = []
    ioFiles          = []

    for f in files:
        if not f.endswith('.csv'):
            notCSVFiles.append(f);
        else:
            if not os.path.exists(f):
                dneFiles.append(f)
            else:
                bname = os.path.basename(f)
                if bname.find('graph') == -1:
                    notGraphFiles.append(f)
                else:
                    outFile = os.path.dirname(f) + '/'
                    outFile += bname.replace('graph', 'cycles')
                    if os.path.exists(outFile):
                        resp = input(outFile + " already exists. Overwrite [y/n]? ")
                        if resp == 'y' or resp == 'Y':
                            ioFiles.append((f, outFile))
                        else:
                            noOverwriteFiles.append(f)
                    else:
                        ioFiles.append((f, outFile))

    print('+-----------------------------------------------------------+')
    print('|                          Results                          |')
    print('+-----------------------------------------------------------+')

    # Initialize lock for printing in parallel
    lock = mp.Lock()
    # Initialize pool for launching tasks asynchronously
    pool = mp.Pool(initializer = init,
                   initargs = (lock,),
                   processes = mp.cpu_count())
    pool.map_async(processFile, ioFiles)
    pool.close()
    pool.join()

    print('+-----------------------------------------------------------+')
    print('|                          Skipped                          |')
    print('+-----------------------------------------------------------+')
    # Print the skipped files nicely


def processFile(ioPair):
    inFile, outFile = ioPair
    cmd = ' '.join([EXE, inFile, outFile])
    subprocess.check_call(shlex.split(cmd))
    lock.acquire()
    print(inFile + ' -> ' + outFile)
    lock.release()


if __name__ == '__main__':
    main(sys.argv[1:]);
