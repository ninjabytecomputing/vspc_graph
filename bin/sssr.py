#!/usr/bin/python3

import os
import itertools
import math
import multiprocessing as mp
import shlex
import subprocess
import sys

# These global variables are changed later
EXE = ''
MAX_LEN = 0

# These global variables are not changed
DONE_STR = 'DONE'
STATUS_STR = 'Status'
OUT_FILE_STR = 'Output Files'
SKIP_FILE_STR = 'Skipped Files'



class Colors:
    End          = '\033[0m'
    Default      = "\033[39m"
    Black        = "\033[30m"
    Red          = "\033[31m"
    Green        = "\033[32m"
    Yellow       = "\033[33m"
    Blue         = "\033[34m"
    Magenta      = "\033[35m"
    Cyan         = "\033[36m"
    LightGray    = "\033[37m"
    DarkGray     = "\033[90m"
    LightRed     = "\033[91m"
    LightGreen   = "\033[92m"
    LightYellow  = "\033[93m"
    LightBlue    = "\033[94m"
    LightMagenta = "\033[95m"
    LightCyan    = "\033[96m"
    White        = "\033[97m"


def init(l):
    global lock
    lock = l


def main(files):
    # Initialize path to the executable, which should be in the same directory
    # as this script if it was installed via CMake
    global EXE
    global MAX_LEN
    EXE = os.path.dirname(os.path.realpath(__file__)) + '/main'

    # Filter out any files with 'cycles' as a substring in the base filename
    files = [f for f in files if 'cycles' not in os.path.basename(f)]

    # Find maximum file length for pretty output
    for f in files:
        if len(f) > MAX_LEN:
            MAX_LEN = len(f)

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

    # Pretty output
    printSeparator()
    printOutputTitle()
    printSeparator()

    # Initialize lock for printing in parallel
    lock = mp.Lock()
    # Initialize pool for launching tasks asynchronously
    pool = mp.Pool(initializer = init,
                   initargs = (lock,),
                   processes = mp.cpu_count())
    # Launch asynchronous processes
    pool.map_async(processFile, ioFiles)
    pool.close()
    pool.join()

    printSeparator()

    print('+-----------------------------------------------------------+')
    print('|                          Skipped                          |')
    print('+-----------------------------------------------------------+')
    # Print the skipped files nicely


def processFile(ioPair):
    inFile, outFile = ioPair
    cmd = ' '.join([EXE, inFile, outFile])
    subprocess.check_call(shlex.split(cmd))
    lock.acquire()
    out = '| ' + inFile
    out += ' ' * (MAX_LEN - len(inFile) + 1)
    out += '|  ' + Colors.Green + DONE_STR + Colors.End + '  |'
    print(out)
    lock.release()


def printSeparator():
    sep  = '+'
    sep += '-' * (MAX_LEN + 2)
    sep += '+'
    sep += '-' * (len(STATUS_STR) + 2)
    sep += '+'
    print(sep)


def printOutputTitle():
    outputTitle = '|'
    outputTitle += ' ' * math.floor((MAX_LEN + 2 - len(OUT_FILE_STR)) / 2)
    outputTitle += OUT_FILE_STR
    outputTitle += ' ' * math.ceil((MAX_LEN + 2 - len(OUT_FILE_STR)) / 2)
    outputTitle += '| ' + STATUS_STR + ' |'
    print(outputTitle)


if __name__ == '__main__':
    main(sys.argv[1:]);
