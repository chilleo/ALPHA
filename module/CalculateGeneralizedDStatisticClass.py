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


class CalculateGeneralizedDStatisticClass(QtCore.QThread):
    def __init__(self, parent=None):
        super(CalculateGeneralizedDStatisticClass, self).__init__(parent)

    def calculate_generalized(self, alignments, species_tree=None, reticulations=None, outgroup=None, window_size=100000000000,
                          window_offset=100000000000, verbose=False, alpha=0.01, use_inv=False, useDir=False,
                 directory="", statistic=False, save=False,  f="DGenStatistic_", plot=False, meta=False, useAlreadyGeneratedStat=True):

        self.emit(QtCore.SIGNAL('GEN_D_10'))

        if(useAlreadyGeneratedStat == False): #generate a dgen stat
            # run the java jar lines goes here

            # FOR REFERENCE, THE ARGS AND JAVA COMMAND
            # String treeString = args[0];
            # String networkString = args[1];
            # String outGroupName = args[2];
            # String saveStatHere = args[3];
            # int numberOfRandomRuns = Integer.parseInt(args[4]);
            # GenerateDgenStatistic(treeString, networkString, outGroupName, saveStatHere, numberOfRandomRuns);

            # Get the global path name to the jar file
            dir_path = os.path.dirname(os.path.realpath(__file__))
            jarPath = os.path.join(dir_path, "DGEN2.jar")

            # Run PhyloNet dgen maker jar file
            numberRandomRuns = 100
            networkString = Create_Network_Helper(species_tree, reticulations, 0.9)
            #species tree and network string need 's to work properly
            species_tree = "'"+species_tree+"'"
            networkString = "'"+networkString+"'"
            jarRunOutput = subprocess.Popen("java -jar {0} {1} {2} {3} {4} {5}".format(jarPath, species_tree, networkString, outgroup, statistic, numberRandomRuns), stdout=subprocess.PIPE,
                                 shell=True)

            # Read output and convert to float
            #pgtst = float(p.stdout.readline())


        self.emit(QtCore.SIGNAL('GEN_D_50'))

        #and then always run the statistic on data (just doing it this way for now to save time. could be chagned later to be slightly more user friendly. but also all users should want to analyze data probably)
        #run the dstat. making temp variables just to keep clear on what is named what (cuz i am changing some things around without messing with some of the code rewriting)
        runInVerboseMode = use_inv
        saveResultsHere = f
        resultsString = run_saved_dgen(statistic, alignments, window_size=window_size, window_offset=window_offset, verbose=runInVerboseMode, alpha=alpha)

        #line here to save results string to file saveResultsHere (do i want to do this or just output to screen?
        # If users want to save the statistic and speed up future runs
        if len(saveResultsHere) > 0:
            num = 0
            file_name = saveResultsHere + ".txt"
            while os.path.exists(file_name):
                file_name = "DGenResults_{0}.txt".format(num)
                num += 1

            with open(file_name, "w") as text_file:
                #output_str = "Taxa: {0}\n".format(taxa)
                #text_file.write(output_str)
                #output_str = "Statistic: {0}\n".format(generate_statistic_string((increase_resized, decrease_resized)))
                #text_file.write(output_str)
                text_file.write(resultsString)
                text_file.close()

        #put a line to either print results or to save em to a file. printing to screen done here
        #self.emit(QtCore.SIGNAL('GEN_D_COMPLETE'))
        self.emit(QtCore.SIGNAL('GEN_D_100'))
        self.emit(QtCore.SIGNAL('DGEN2_FINISHED'), resultsString)

        debugHere = 0

        #run_saved_dgen(?,
        #               ['/Users/leo/rice/res/data/cichlid/alignment/cichlid6tax.phylip-sequential.txt'],
        #               verbose=True, plot='/Users/leo/rice/res/data/dgen/tmp/figC/plot_figCVerbose', meta='Dgen')

        # OLD WAY
        # alignments_to_d_resized, alignments_to_windows_to_d, standard_o, verbose_o = CalculateGeneralizedDStatistic.calculate_generalized\
        #     (alignments, species_tree, reticulations, outgroup, window_size, window_offset, verbose, alpha, use_inv,
        #                                                     useDir, directory, statistic, save, f, plot, meta)

        # self.emit(QtCore.SIGNAL("L_FINISHED"), alignments_to_d_resized, alignments_to_windows_to_d, standard_o, verbose_o)

    def run(self):
        """
            Starts PyQt Thread. Called with "start()".
        """
        # try:
        #     self.window_splitter(self.inputFilename, self.windowSize, self.windowOffset)
        # except IOError:
        #     self.emit(QtCore.SIGNAL('INVALID_ALIGNMENT_FILE'), self.inputFilename)
        #     return

        self.calculate_generalized(self.alignments,
                                   species_tree=self.species_tree,
                                   reticulations=self.r,
                                   outgroup=self.o,
                                   window_size=self.window_size,
                                   window_offset=self.window_offset,
                                   verbose=self.verbose,
                                   alpha=self.alpha,
                                   use_inv=self.use_inv,
                                   useDir=self.useDir,
                                   directory=self.directory,
                                   statistic=self.statistic,
                                   save=self.save,
                                   f=self.save_location,
                                   plot=self.plot,
                                   meta=self.meta,
                                   useAlreadyGeneratedStat=self.useAlreadyGeneratedStat)

        #self.emit(QtCore.SIGNAL('GEN_D_COMPLETE'), None)

if __name__ == '__main__':
    gd = CalculateGeneralizedDStatisticClass()

    species_tree = '((P1,P2),(P3,O));'
    # species_tree = '(((P1,P2),(P3,(P4,P5))),O);'
    r = [('P3', 'P1')]
    alignments = ["exampleFiles/seqfile.txt"]

    if platform == "darwin":
        alignments = ["/Users/Peter/PycharmProjects/ALPHA/exampleFiles/seqfile.txt"]
    else:
        alignments = ["C:\\Users\\travi\Desktop\\dFoilStdPlusOneFar50kbp\\dFoilStdPlusOneFar50kbp\\sim2\\seqfile.txt"]

    # gd.calculate_generalized(alignments, species_tree, r, window_size=50000, window_offset=50000, verbose=True, alpha=0.01, save=True)
    gd.calculate(alignments, species_tree, r, outgroup="O", window_size=50000, window_offset=50000, verbose=True, alpha=0.01, save=True)

    # save_file = "C:\\Users\\travi\\Documents\\ALPHA\\CommandLineFiles\\DGenStatistic_35.txt"
    # plot_formatting(calculate_generalized(alignments, statistic=save_file))

    # print calculate_generalized(alignments, statistic="C:\\Users\\travi\\Documents\\ALPHA\\CommandLineFiles\\DGenStatistic_10.txt", verbose=True)
    # calculate_generalized(alignments, statistic="C:\\Users\\travi\\Documents\\ALPHA\\CommandLineFiles\\DGenStatistic_35.txt")

    # python - c "from CalculateGeneralizedDStatistic import *; calculate_generalized(['C:\\Users\\travi\\Documents\\PhyloVis\\exampleFiles\\ExampleDFOIL.phylip'], statistic='C:\\Users\\travi\\Documents\\ALPHA\\CommandLineFiles\\DGenStatistic_35.txt')"

    # species_tree, r = '(((P1,P2),(P3,(P4,P5))),O);', [('P1', 'P3')]
    # alignments = ["C:\\Users\\travi\\Documents\\PhyloVis\\exampleFiles\\ExampleDFOIL.phylip"]
    # alignments = ["C:\\Users\\travi\\Desktop\\sixtaxa.txt"]
    # i = calculate_generalized(alignments, species_tree, r, 100000, 100000, True, save=True)

    # for j in range(10):
    #     k = calculate_generalized(alignments, species_tree, r, 100000, 100000, True, save=True)
    #     if i != k:
    #         print "FAIL"
    #         print i
    #         print k
    #     print j


    # print pattern_string_generator(['A', 'A', 'A', 'A', 'A'])

    # Inputs for paper
    # file = "C:\\Users\\travi\\Desktop\\concatFile.phylip.txt"
    # species_tree = '((C,G),(((A,Q),L),R));'
    #
    # window_size, window_offset = 10000, 1000
    # r = [('L', 'R')]
    # plot_formatting(calculate_generalized(file, species_tree, r, window_size, window_offset, True))
    # window_size, window_offset = 100000, 10000
    # plot_formatting(calculate_generalized(file, species_tree, r, window_size, window_offset, True))
    #
    # window_size, window_offset = 10000, 1000
    # r = [('Q', 'R')]
    # plot_formatting(calculate_generalized(file, species_tree, r, window_size, window_offset, True))
    # window_size, window_offset = 100000, 10000
    # plot_formatting(calculate_generalized(file, species_tree, r, window_size, window_offset, True))
    #
    # window_size, window_offset = 10000, 1000
    # r = [('Q', 'G')]
    # plot_formatting(calculate_generalized(file, species_tree, r, window_size, window_offset, True))
    # window_size, window_offset = 100000, 10000
    # plot_formatting(calculate_generalized(file, species_tree, r, window_size, window_offset, True))

    # concat_directory("/Users/Peter/PycharmProjects/ALPHA/test_phylip_dir")
    # print calculate_generalized('/Users/Peter/PycharmProjects/ALPHA/CLFILE', '(((P1,P2),(P3,P4)),O);', [('P1', 'P3')], 50000, 50000, True)


    # file = 'C:\\Users\\travi\\Desktop\\clphylipseq.txt'
    # # r = [('L', 'R')]
    # r = [('Q', 'R')]
    # # r = [('Q', 'G')]
    # print calculate_generalized(file , '((C,G),(((A,Q),L),R));', r, 100000, 100000, True)

    # concat_directory("/Users/Peter/PycharmProjects/ALPHA/travy_test")
    # print calculate_generalized('/Users/Peter/PycharmProjects/ALPHA/CLFILE', '(((P1,P2),(P3,P4)),O);', [('P1', 'P3')], 50000, 50000, True)

    # plot_formatting(calculate_generalized(alignments, species_tree, r, 1000, 1000, True))
    # # lstat, signif, windows_to_l = calculate_generalized(alignment, species_tree, r, 1000, 1000, True, 0.05)
    # # plot_formatting((lstat, signif, windows_to_l))
    # plot_formatting(calculate_generalized('C:\\Users\\travi\\Desktop\\seqfileNamed', '(((P1,P2),(P3,P4)),O);', [('P3', 'P1')], 1000, 1000, False, 0.99), False)

    # print calculate_generalized('C:\\Users\\travi\\Desktop\\seqfileNamed', '(((P1,P2),(P3,P4)),O);', [('P1', 'P3')], 50000, 50000, True)

# python -c"from CalculateGeneralizedDStatistic import *; plot_formatting(calculate_generalized('C:\\Users\\travi\\Desktop\\seqfileNamed', '(((P1,P2),(P3,P4)),O);', [('P1', 'P3')], 100000, 100000, True, 0.01), True)"
