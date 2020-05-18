#!/usr/bin/python3

import os
import shlex
import subprocess
import sys


EXECUTABLE = './release/src/main'

def main(files):
    # Filter out any files with '_cycles' as a substring in the filename
    files = [f for f in files if '_cycles' not in f]

    skippedFiles = [];
    fileMap = {};

    for f in files:
        if not f.endswith('csv'):
            skippedFiles.append(f);
        else:
            if not os.path.exists(f):
                print(f + " doesn't exist - skipping")
            else:
                outFile = f[:-4] + '_cycles.csv'
                if os.path.exists(outFile):
                    resp = input(outFile + " already exists. Overwrite [y/n]? ")
                    if resp == 'y' or resp == 'Y':
                        fileMap[f] = outFile
                else:
                    fileMap[f] = outFile

    print('+-----------------------------------------------------------+')
    print('|                          Results                          |')
    print('+-----------------------------------------------------------+')
    for inF, outF in fileMap.items():
        cmd = ' '.join([EXECUTABLE, inF, outF])
        subprocess.check_call(shlex.split(cmd))


if __name__ == '__main__':
    main(sys.argv[1:]);
