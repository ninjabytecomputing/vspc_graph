#!/usr/bin/python3

import os
import shlex
import subprocess
import sys

EXE = './release/src/main'


def main(files):
    # Filter out any files with 'cycle' as a substring in the base filename
    files = [f for f in files if 'cycle' not in os.path.basename(f)]

    notCSVFiles      = []
    dneFiles         = []
    notGraphFiles    = []
    noOverwriteFiles = []
    fileMap          = {}

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
                    outFile += bname.replace('graph', 'cycle')
                    if os.path.exists(outFile):
                        resp = input(outFile + " already exists. Overwrite [y/n]? ")
                        if resp == 'y' or resp == 'Y':
                            fileMap[f] = outFile
                        else:
                            noOverwriteFiles.append(f)
                    else:
                        fileMap[f] = outFile

    print('+-----------------------------------------------------------+')
    print('|                          Results                          |')
    print('+-----------------------------------------------------------+')
    for inF, outF in fileMap.items():
        cmd = ' '.join([EXE, inF, outF])
        # subprocess.check_call(shlex.split(cmd))

    print('+-----------------------------------------------------------+')
    print('|                          Skipped                          |')
    print('+-----------------------------------------------------------+')
    # Print the skipped files nicely


if __name__ == '__main__':
    main(sys.argv[1:]);
