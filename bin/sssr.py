#!/usr/bin/python3

import math
import multiprocessing as mp
import os
import shlex
import subprocess
import sys

# This will be initialized during main.
EXE = ''


class FileManager:
    CONST_OUT_FILE  = 'Output Files'
    CONST_SKIP_FILE = 'Skipped Files'
    CONST_DONE      = 'Done'
    CONST_STATUS    = 'Status'
    CONST_DNE       = 'DNE'
    CONST_NO_WRITE  = 'No write'

    def __init__(self, files):
        # Keep csv files
        files = [f for f in files if os.path.basename(f).endswith('csv')]
        # Keep files whose base file name does not begin with 'cycles'
        files = [f for f in files if not os.path.basename(f).startswith('cycles')]
        # Keep files whose base file name starts with graph
        files = [f for f in files if os.path.basename(f).startswith('graph')]

        self.filesDNE         = []
        self.filesNoOverwrite = []
        self.filesIO          = []

        self.maxFileLength = max(len(self.CONST_OUT_FILE),
                                 len(self.CONST_SKIP_FILE))

        self.separator    = ''
        self.outputTitle  = ''
        self.skippedTitle = ''

        for f in files:
            self._sortFile(f)

        self._initSeparator()
        self._initOutputTitle()
        self._initSkippedTitle()

    def launchParallel(self):
        # Early out if there are no files to process
        if len(self.filesIO) > 0:
            self._printSeparator()
            self._printOutputTitle()
            self._printSeparator()

            # Initialize lock for printing in parallel
            lock = mp.Lock()
            # Initialize pool for launching tasks asynchronously
            pool = mp.Pool(initializer = self._initLock,
                           initargs    = (lock,),
                           processes   = mp.cpu_count())
            # Launch asynchronous processes
            pool.map_async(self._processFile, self.filesIO)
            pool.close()
            pool.join()
            self._printSeparator()

        if len(self.filesDNE) > 0 or len(self.filesNoOverwrite) > 0:
            # If there were no files to process, then we need to first
            # print a separator line.
            if len(self.filesIO) == 0:
                self._printSeparator()
            self._printSkippedTitle()
            self._printSeparator()

            # Print skipped files
            self._printSkippedFiles()
            self._printSeparator()


    def _printSeparator(self):
        print(self.separator)

    def _printOutputTitle(self):
        print(self.outputTitle)

    def _printSkippedTitle(self):
        print(self.skippedTitle)

    def _printSkippedFiles(self):
        # Red status for files that don't exist
        for f in self.filesDNE:
            out = '| ' + f
            out += ' ' * (self.maxFileLength - len(f) + 1)
            out += '|   ' + '\033[91m' + self.CONST_DNE + '\033[0m' + '    |'
            print(out)
        # Cyan status for files marked as no overwrite
        for f in self.filesNoOverwrite:
            out = '| ' + f
            out += ' ' * (self.maxFileLength - len(f) + 1)
            out += '| ' + '\033[96m' + self.CONST_NO_WRITE + '\033[0m' + ' |'
            print(out)

    def _processFile(self, ioPair):
        greenStart = '\033[32m'
        greenEnd   = '\033[0m'
        inFile, outFile = ioPair
        cmd = ' '.join([EXE, inFile, outFile])
        subprocess.check_call(shlex.split(cmd))
        lock.acquire()
        stat = '| ' + outFile
        stat += ' ' * (self.maxFileLength - len(outFile) + 1)
        stat += '|   ' + '\033[32m' + self.CONST_DONE + '\033[0m' + '   |'
        print(stat)
        lock.release()

    def _sortFile(self, f):
        if not os.path.exists(f):
            if len(f) > self.maxFileLength:
                self.maxFileLength = len(f)
            self.filesDNE.append(f)
        else:
            outFile = os.path.dirname(f) + '/'
            outFile += os.path.basename(f).replace('graph', 'cycles')
            if os.path.exists(outFile):
                resp = input(outFile + " already exists. Overwrite [y/n]? ")
                if resp == 'y' or resp == 'Y':
                    if len(outFile) > self.maxFileLength:
                        self.maxFileLength = len(outFile)
                    self.filesIO.append((f, outFile))
                else:
                    if len(f) > self.maxFileLength:
                        self.maxFileLength = len(f)
                    self.filesNoOverwrite.append(f)
            else:
                if len(outFile) > self.maxFileLength:
                    self.maxFileLength = len(outFile)
                self.filesIO.append((f, outFile))

    def _initLock(self, l):
        global lock
        lock = l

    def _initSeparator(self):
        # Note: CONST_NOT_GRAPH is the longest status name
        self.separator  = '+'
        self.separator += '-' * (self.maxFileLength + 2)
        self.separator += '+'
        self.separator += '-' * (len(self.CONST_NO_WRITE) + 2)
        self.separator += '+'

    def _initOutputTitle(self):
        totalSpaceLen = self.maxFileLength + 2 - len(self.CONST_OUT_FILE)
        preSpaceLen   = math.floor((totalSpaceLen) / 2)
        postSpaceLen  = totalSpaceLen - preSpaceLen
        self.outputTitle = '|'
        self.outputTitle += ' ' * preSpaceLen
        self.outputTitle += self.CONST_OUT_FILE
        self.outputTitle += ' ' * postSpaceLen
        self.outputTitle += '|  ' + self.CONST_STATUS + '  |'

    def _initSkippedTitle(self):
        totalSpaceLen = self.maxFileLength + 2 - len(self.CONST_SKIP_FILE)
        preSpaceLen   = math.floor((totalSpaceLen) / 2)
        postSpaceLen  = totalSpaceLen - preSpaceLen
        self.skippedTitle = '|'
        self.skippedTitle += ' ' * preSpaceLen
        self.skippedTitle += self.CONST_SKIP_FILE
        self.skippedTitle += ' ' * postSpaceLen
        self.skippedTitle += '|  ' + self.CONST_STATUS + '  |'


def main(files):
    # Initialize path to the executable, which should be in the same directory
    # as this script if it was installed via CMake
    global EXE
    EXE = os.path.dirname(os.path.realpath(__file__)) + '/main'

    # Initialize FileManager class, which does everything for us
    manager = FileManager(files)
    # Launch processes in parallel
    manager.launchParallel()


if __name__ == '__main__':
    main(sys.argv[1:]);
