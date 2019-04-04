from PyQt4 import QtCore

from sys import *
import os
import subprocess

# if platform == 'win32':
#     path.insert(0, "../CommandLineFiles")
# elif platform == 'darwin':
#     path.insert(0, "..\\CommandLineFiles")

# import CalculateGeneralizedDStatistic
import sys

from CommandLineFiles.RunDGEN import run_saved_dgen, Create_Network_Helper

sys.path.append('..\\')

from CommandLineFiles import CalculateGeneralizedDStatistic

"""
Functions:

~
Chabrielle Allen
Travis Benedict
Peter Dulworth
"""


class RunSW(QtCore.QThread):
    def __init__(self, parent=None):
        super(RunSW, self).__init__(parent)

    def run_sw_jar(self, sequence_path, sequence_length, window_size=100000000000,
                          window_offset=100000000000):


        self.emit(QtCore.SIGNAL('SW_UPDATE'),
                  "Running Smooth Winds....")

        # run the java jar lines goes here

        # FOR REFERENCE, THE ARGS AND JAVA COMMAND
        # String seq path = args[0];
        # String seq length = args[1];
        # String window size = args[2];
        # String window offset = args[3];

        # Get the global path name to the jar file
        dir_path = os.path.dirname(os.path.realpath(__file__))
        jarPath = os.path.join(dir_path, "SmoothWinds.jar")

        # Run SW jar file

        #UNCOMMENT ME TO RUN JAR
        jarRunOutput = subprocess.Popen("java -jar {0} {1} {2} {3} {4}".format(jarPath, sequence_path, sequence_length, window_size, window_offset), stdout=subprocess.PIPE,
                             shell=True)

        # Read output and convert to float
        #pgtst = float(p.stdout.readline())
        #waitTilJarIsDoneRunning = p.stdout.readline()
        jarRunOutput.wait()


        self.emit(QtCore.SIGNAL('SW_UPDATE'), "Run Complete. Output stored in SmoothWinds folder created in the same directory as the fasta file.")

        debugHere = 0


    def run(self):
        """
            Starts PyQt Thread. Called with "start()".
        """
        # try:
        #     self.window_splitter(self.inputFilename, self.windowSize, self.windowOffset)
        # except IOError:
        #     self.emit(QtCore.SIGNAL('INVALID_ALIGNMENT_FILE'), self.inputFilename)
        #     return

        self.run_sw_jar(sequence_path=self.sequencePathText, sequence_length=self.sequenceLengthFloat, window_size=self.windowSizeFloat,window_offset=self.windowOffsetFloat)

        #self.run_sw_jar(sequence_path=self.sequence_path,
        #                   sequence_length=self.sequence_length,
        #                   window_size=self.window_size,
        #                   window_offset=self.window_offset)

        #self.emit(QtCore.SIGNAL('GEN_D_COMPLETE'), None)
