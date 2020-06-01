#!/usr/bin/python3

import argparse
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

    def __init__(self, files, overwrite, verbose):
        # Keep csv files
        files = [f for f in files if os.path.basename(f).endswith('csv')]
        # Keep files whose base file name does not begin with 'cycles'
        files = [f for f in files if not os.path.basename(f).startswith('cycles')]
        # Keep files whose base file name starts with graph
        files = [f for f in files if os.path.basename(f).startswith('graph')]


        self.filesDNE         = []
        self.filesNoOverwrite = []
        self.filesIO          = []

        self.verbose       = verbose
        self.maxFileLength = max(len(self.CONST_OUT_FILE),
                                 len(self.CONST_SKIP_FILE))

        self.separator    = ''
        self.outputTitle  = ''
        self.skippedTitle = ''

        for f in files:
            self._sortFile(f, overwrite)

        self._initSeparator()
        self._initOutputTitle()
        self._initSkippedTitle()

    def launchParallel(self):
        '''
        Process all files in parallel. If verbose output is enabled,
        status updates are printed to the terminal when each file has
        finished. After all files have been processed, all of the
        input files that were skipped are then printed as well.
        '''
        # Early out if there are no files to process
        if len(self.filesIO) > 0:
            if self.verbose:
                self._printSeparator()
                self._printOutputTitle()
                self._printSeparator()

            # Initialize lock for printing in parallel
            lock = mp.Lock()
            # Initialize pool for launching tasks asynchronously
            pool = mp.Pool(initializer = self._initLock,
                           initargs    = (lock,),
                           processes   = mp.cpu_count()-1)
            # Launch asynchronous processes
            pool.map_async(self._processFile, self.filesIO)
            pool.close()
            pool.join()
            if self.verbose:
                self._printSeparator()

        if self.verbose:
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
        ''' Print separator '''
        print(self.separator)

    def _printOutputTitle(self):
        ''' Print title for the output files section '''
        print(self.outputTitle)

    def _printSkippedTitle(self):
        ''' Print title for the skipped files section '''
        print(self.skippedTitle)

    def _printSkippedFiles(self):
        ''' Print all skipped files. '''
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
        '''
        Call the C++ executable on the input-output file name pair, and print
        a status update when it's done.
        '''
        greenStart = '\033[32m'
        greenEnd   = '\033[0m'
        inFile, outFile = ioPair
        cmd = ' '.join([EXE, inFile, outFile])
        subprocess.check_call(shlex.split(cmd))
        if self.verbose:
            lock.acquire()
            stat = '| ' + outFile
            stat += ' ' * (self.maxFileLength - len(outFile) + 1)
            stat += '|   ' + '\033[32m' + self.CONST_DONE + '\033[0m' + '   |'
            print(stat)
            lock.release()

    def _sortFile(self, f, overwrite):
        '''
        Sort the input file into one of three categories for processing:
            1) input file doesn't exist,
            2) output file exists but shouldn't be overwritten, or
            3) input file should be processed.
        '''
        if not os.path.exists(f):
            if len(f) > self.maxFileLength:
                self.maxFileLength = len(f)
            self.filesDNE.append(f)
        else:
            outFile = os.path.dirname(f) + '/'
            outFile += os.path.basename(f).replace('graph', 'cycles')
            if os.path.exists(outFile) and not overwrite:
                if len(f) > self.maxFileLength:
                    self.maxFileLength = len(f)
                self.filesNoOverwrite.append(f)
            else:
                if len(outFile) > self.maxFileLength:
                    self.maxFileLength = len(outFile)
                self.filesIO.append((f, outFile))

    def _initLock(self, l):
        '''
        Initialize the global lock that is used to write the output in
        a thread-safe manner.
        '''
        global lock
        lock = l

    def _initSeparator(self):
        '''
        Construct the string for the horizontal bar (separator) for the
        output.
        '''
        # Note: CONST_NOT_GRAPH is the longest status name
        self.separator  = '+'
        self.separator += '-' * (self.maxFileLength + 2)
        self.separator += '+'
        self.separator += '-' * (len(self.CONST_NO_WRITE) + 2)
        self.separator += '+'

    def _initOutputTitle(self):
        '''
        Construct the string for the title of the 'output files' section
        of the output.
        '''
        totalSpaceLen = self.maxFileLength + 2 - len(self.CONST_OUT_FILE)
        preSpaceLen   = math.floor((totalSpaceLen) / 2)
        postSpaceLen  = totalSpaceLen - preSpaceLen
        self.outputTitle = '|'
        self.outputTitle += ' ' * preSpaceLen
        self.outputTitle += self.CONST_OUT_FILE
        self.outputTitle += ' ' * postSpaceLen
        self.outputTitle += '|  ' + self.CONST_STATUS + '  |'

    def _initSkippedTitle(self):
        '''
        Construct the string for the title of the 'skipped files' section
        of the output.
        '''
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

    # Check for any flags
    parser = argparse.ArgumentParser(
        description = 'Process graph connectivity data to detect all cycles. \
            The connectivity data must be supplied in a CSV file where each \
            row represents an undirected edge by listing the indices of the \
            two nodes that flank the edge.')
    parser.add_argument('files', type=str, nargs='+',
                        help='Files to be processed')
    parser.add_argument('--overwrite', action='store_true',
                        help='overwrite all output files already generated')
    parser.add_argument('--verbose', action='store_true',
                        help='print file information during runtime')
    args = parser.parse_args()

    # Initialize FileManager class, which does everything for us
    manager = FileManager(args.files, args.overwrite, args.verbose)
    # Launch processes in parallel
    manager.launchParallel()


if __name__ == '__main__':
    main(sys.argv[1:]);
