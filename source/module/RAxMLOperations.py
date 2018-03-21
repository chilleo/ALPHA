from natsort import natsorted
from sys import platform
import subprocess
import shutil
import os
import re
from PyQt4 import QtCore
from Bio import Phylo

"""
Functions:
    __init__(self, parent=None)
    taxon_names_getter(self, phylip)
    raxml_species_tree(self, phylip, rooted=False, outgroup=None, customRax=False, customRaxCommand='', output_directory="RAxML_SpeciesTree")
    rooter(self, newick_file, outgroup)
    window_splitter(self, filename, window_size, step_size)
    raxml_windows(self, numBootstraps, model, rooted=False, outgroup=None, window_directory='windows', output_directory='RAxML_Files')
    run(self)
~
Chabrielle Allen
Travis Benedict
Peter Dulworth
"""


class RAxMLOperations(QtCore.QThread):
    def __init__(self, parent=None):
        super(RAxMLOperations, self).__init__(parent)

    def taxon_names_getter(self, phylip):
        """
        Creates a list of taxon names from the inputted file
        Input:
        phylip --- a file inputted by the user
        Output:
        taxon_names --- a list of taxon names from the inputted phylip file
        """

        # Initialize a list for the taxa
        taxon_names = []

        with open(phylip) as f:
            # Create a list of each line in the file
            lines = f.readlines()

        # Iterate over each line after the first one
        for line in lines[1:]:
            # Add each sequence to a list
            taxon = line.split()[0]
            taxon_names.append(taxon)

        return taxon_names

    def raxml_species_tree(self, phylip, rooted=False, outgroup=None, customRax=False, customRaxCommand='', output_directory="RAxML_SpeciesTree"):
        """
        Runs RAxML on input PHYLIP file to create a species
        tree.

        Inputs:
        phylip -- a file inputted by the user.

        Returns:
        A species tree folder.
        """

        # Delete the folder and remake it if it already exists
        if os.path.exists(output_directory):
            shutil.rmtree(output_directory)
        os.makedirs(output_directory)
        self.emit(QtCore.SIGNAL('SPECIES_TREE_PER'), 10)

        if customRax:
            # run RAxML & wait until command line is finished running
            p = subprocess.Popen(customRaxCommand + " -s {0} -n txt".format(phylip), shell=True)
            p.wait()
        else:
            # run RAxML & wait until command line is finished running
            p = subprocess.Popen("raxmlHPC -f a -x12345 -p 12345 -# 2 -m GTRGAMMA -s {0} -n txt".format(phylip), shell=True)
            p.wait()

        self.emit(QtCore.SIGNAL('SPECIES_TREE_PER'), 60)

        # Regular expression for identifying floats
        float_pattern = "([+-]?\\d*\\.\\d+)(?![-+0-9\\.])"

        # Create a separate file with the topology of the best tree
        with open("RAxML_bestTree.txt") as f:
            # Read newick string from file
            topology = f.readline()
            # Delete float branch lengths, ":" and "\n" from newick string
            topology = ((re.sub(float_pattern, '', topology)).replace(":", "")).replace("\n", "")
            file = open("Topology_bestTree.txt", "w")
            file.write(topology)
            file.close()
        self.emit(QtCore.SIGNAL('SPECIES_TREE_PER'), 80)

        # If rooting is desired root the appropriate files
        if rooted:
            self.rooter("RAxML_bestTree.txt", outgroup)
            # self.rooter("RAxML_result.txt", outgroup)

        # windows
        if platform == "win32":
            # Move RAxML output files into their own destination folder - Windows
            os.rename("RAxML_bestTree.txt", output_directory + "\RAxML_ST_bestTree.txt")
            os.rename("RAxML_bipartitions.txt", output_directory + "\RAxML_ST_bipartitions.txt")
            os.rename("RAxML_bipartitionsBranchLabels.txt", output_directory + "\RAxML_ST_bipartitionsBranchLabels.txt")
            os.rename("RAxML_bootstrap.txt", output_directory + "\RAxML_ST_bootstrap.txt")
            os.rename("RAxML_info.txt", output_directory + "\RAxML_ST_info.txt")
            os.rename("topology_bestTree.txt", output_directory + "\Topology_ST_bestTree.txt")

        # mac
        elif platform == "darwin":
            # Move RAxML output files into their own destination folder - Mac
            os.rename("RAxML_bestTree.txt", output_directory + "/RAxML_ST_bestTree.txt")
            os.rename("RAxML_bipartitions.txt", output_directory + "/RAxML_ST_bipartitions.txt")
            os.rename("RAxML_bipartitionsBranchLabels.txt", output_directory + "/RAxML_ST_bipartitionsBranchLabels.txt")
            os.rename("RAxML_bootstrap.txt", output_directory + "/RAxML_ST_bootstrap.txt")
            os.rename("RAxML_info.txt", output_directory + "/RAxML_ST_info.txt")
            os.rename("topology_bestTree.txt", output_directory + "/Topology_ST_bestTree.txt")

        if platform == 'win32':
            with open('RAxML_SpeciesTree\RAxML_ST_bestTree.txt', 'r') as f:
                self.speciesTree = f.read().replace('\n', '')
        elif platform == 'darwin':
            with open('RAxML_SpeciesTree/RAxML_ST_bestTree.txt', 'r') as f:
                self.speciesTree = f.read().replace('\n', '')

        self.emit(QtCore.SIGNAL('SPECIES_TREE_PER'), 100)
        self.emit(QtCore.SIGNAL('SPECIES_TREE_COMPLETE'), 'Species Tree Generated', "'Show Details...' to view the species tree newick.", self.speciesTree)
        self.emit(QtCore.SIGNAL('SPECIES_TREE_COMPLETE_RETURN_ST'), self.speciesTree)
        print self.speciesTree

    def rooter(self, newick_file, outgroup):
        """
        Rewrites tree newick strings to be rooted
        Inputs:
        newick_file --- the file containing the newick string
        outgroup --- the outgroup to root at
        """

        # Create the tree object and root it
        tree = Phylo.read(newick_file, "newick")
        tree.rooted = True
        tree.root_with_outgroup(outgroup)

        # Write the newick string over the previous one
        Phylo.write(tree, newick_file, "newick")

    def window_splitter(self, filename, window_size, step_size):
        """
        Creates smaller PHYLIP files based on a window size inputted into
        the GUI.

        Inputs:
        filename --- name of the PHYLIP file to be used
        window_size --- the number of nucleotides to include in each window
        step_size --- the number of nucleotides between the beginning of each window

        Output:
        Smaller "window" files showing sections of the genome in PHYLIP format.
        """

        output_folder = "windows"

        # Delete the folder and remake it
        if os.path.exists(output_folder):
            shutil.rmtree(output_folder)

        # Use os.system instead of os.makedirs(output_folder) so that we can sudo and avoid permissions errors when deploying
        os.system("sudo mkdir " + output_folder)

        # Create a list for the output files
        output_files = []

        with open(filename) as f:
            # First line contains the number and length of the sequences
            line = f.readline()
            line = line.split()

            number_of_sequences = int(line[0])
            length_of_sequences = int(line[1])

            # Initialize a pointer for the beginning of each window
            i = 0
            # Initialize a count for the total number of windows
            num_windows = 0

            # Determine the total number of windows needed
            while (i + window_size - 1 < length_of_sequences):
                i += step_size
                num_windows += 1

            # Create a file for each window and add it to the list
            for i in range(num_windows):
                output_files.append(open(output_folder + "/window" + str(i) + ".phylip", "w"))
                output_files[i].close()

            # Write the number and length of the sequences to each file
            for i in range(num_windows):
                file = open(output_files[i].name, "a")
                file.write(" " + str(number_of_sequences) + " ")
                file.write(str(window_size) + "\n")
                file.close()

            # Subsequent lines contain taxon and sequence separated by a space
            for i in range(number_of_sequences):
                line = f.readline()
                line = line.split()

                taxon = line[0]
                sequence = line[1]

                for j in range(num_windows):
                    l = j * step_size
                    file = open(output_files[j].name, "a")
                    file.write(taxon + " ")
                    window = ""
                    for k in range(window_size):
                        window += sequence[l + k]

                    # Writes file to folder
                    file.write(window + "\n")
                    file.close()

    def raxml_windows(self, numBootstraps, model, rooted=False, outgroup=None, window_directory='windows', output_directory='RAxML_Files'):
        """
        Runs RAxML on files in the directory containing files from
        window_splitter().

        Inputs:
        numBootstraps --- the number of bootstraps to use
        model --- the type of model to use in RAxML
        """

        topology_output_directory = "Topologies"

        # Delete the folder and remake it if it already exists
        if os.path.exists(output_directory):
            shutil.rmtree(output_directory)

        # Use os.system instead of os.makedirs(output_directory) so that we can sudo and avoid permissions errors when deploying
        os.system("sudo mkdir " + output_directory)

        # Delete the folder and remake it if it already exists
        if os.path.exists(topology_output_directory):
            shutil.rmtree(topology_output_directory)

        # Use os.system instead of os.makedirs(topology_output_directory) so that we can sudo and avoid permissions errors when deploying
        os.system("sudo mkdir " + topology_output_directory)

        percent_complete = 0

        # Iterate over each folder in the given directory in numerical order
        for filename in natsorted(os.listdir(window_directory)):

            # If file is a phylip file run RAxML on it
            if filename.endswith(".phylip"):

                file_number = filename.replace("window", "")
                file_number = file_number.replace(".phylip", "")

                input_file = os.path.join(window_directory, filename)

                # Run RAxML
                # NOT custom raxml
                if not self.isCustomRaxmlCommand:
                    if self.bootstrap:
                        p = subprocess.Popen("raxmlHPC -f a -x12345 -p 12345 -# {2} -m {3} -s {0} -n {1}".format(input_file, file_number, numBootstraps, model), shell=True)
                    else:
                        p = subprocess.Popen("raxmlHPC -d -p 12345 -m {2} -s {0} -n {1}".format(input_file, file_number, model), shell=True)
                # custom raxml
                else:
                    p = subprocess.Popen(self.customRaxmlCommand + " -s {0} -n {1}".format(input_file, file_number), shell=True)

                # Wait until command line is finished running
                p.wait()

                # Regular expression for identifying floats
                float_pattern = "([+-]?\\d*\\.\\d+)(?![-+0-9\\.])"

                # Create a separate file with the topology of the best tree
                with open("RAxML_bestTree." + file_number) as f:
                    # Read newick string from file
                    topology = f.readline()

                    # Delete float branch lengths, ":" and "\n" from newick string
                    topology = ((re.sub(float_pattern, '', topology)).replace(":", "")).replace("\n", "")

                    # Otherwise write the topology newick string to a file
                    # else:
                    if platform == 'win32':
                        file = open(topology_output_directory + "\Topology_bestTree." + file_number, "w")
                        file.write(topology)
                        file.close()
                    elif platform == 'darwin':
                        file = open(topology_output_directory + "/Topology_bestTree." + file_number, "w")
                        file.write(topology)
                        file.close()

                if self.bootstrap:

                    # If rooting is desired root the appropriate files
                    if rooted:
                        self.rooter("RAxML_bestTree." + file_number, outgroup)
                        self.rooter("RAxML_bipartitions." + file_number, outgroup)
                        self.rooter("RAxML_bipartitionsBranchLabels." + file_number, outgroup)

                    if platform == "win32":
                        # Move RAxML output files into their own destination folder - Windows
                        os.rename("RAxML_bestTree." + file_number, output_directory + "\RAxML_bestTree." + file_number)
                        os.rename("RAxML_bipartitions." + file_number, output_directory + "\RAxML_bipartitions." + file_number)
                        os.rename("RAxML_bipartitionsBranchLabels." + file_number, output_directory + "\RAxML_bipartitionsBranchLabels." + file_number)
                        os.rename("RAxML_bootstrap." + file_number, output_directory + "\RAxML_bootstrap." + file_number)
                        os.rename("RAxML_info." + file_number, output_directory + "\RAxML_info." + file_number)
                        # os.rename("topology_bestTree." + file_number, topology_output_directory + "\Topology_bestTree." + file_number)

                    elif platform == "darwin":
                        # Move RAxML output files into their own destination folder - Mac
                        os.rename("RAxML_bestTree." + file_number, output_directory + "/RAxML_bestTree." + file_number)
                        os.rename("RAxML_bipartitions." + file_number, output_directory + "/RAxML_bipartitions." + file_number)
                        os.rename("RAxML_bipartitionsBranchLabels." + file_number, output_directory + "/RAxML_bipartitionsBranchLabels." + file_number)
                        os.rename("RAxML_bootstrap." + file_number, output_directory + "/RAxML_bootstrap." + file_number)
                        os.rename("RAxML_info." + file_number, output_directory + "/RAxML_info." + file_number)
                        # os.rename("topology_bestTree." + file_number, topology_output_directory + "/Topology_bestTree." + file_number)

                else:

                    # If rooting is desired root the appropriate files
                    if rooted:
                        self.rooter("RAxML_bestTree." + file_number, outgroup)
                        self.rooter("RAxML_result." + file_number, outgroup)

                    if platform == "win32":
                        # Move RAxML output files into their own destination folder - Windows
                        os.rename("RAxML_bestTree." + file_number, output_directory + "\RAxML_bestTree." + file_number)
                        os.rename("RAxML_log." + file_number, output_directory + "\RAxML_log." + file_number)
                        os.rename("RAxML_randomTree." + file_number, output_directory + "\RAxML_randomTree." + file_number)
                        os.rename("RAxML_result." + file_number, output_directory + "\RAxML_result." + file_number)
                        os.rename("RAxML_info." + file_number, output_directory + "\RAxML_info." + file_number)
                        # os.rename("topology_bestTree." + file_number, topology_output_directory + "\Topology_bestTree." + file_number)

                    elif platform == "darwin":
                        # Move RAxML output files into their own destination folder - Mac
                        os.rename("RAxML_bestTree." + file_number, output_directory + "/RAxML_bestTree." + file_number)
                        os.rename("RAxML_log." + file_number, output_directory + "/RAxML_log." + file_number)
                        os.rename("RAxML_randomTree." + file_number, output_directory + "/RAxML_randomTree." + file_number)
                        os.rename("RAxML_result." + file_number, output_directory + "/RAxML_result." + file_number)
                        os.rename("RAxML_info." + file_number, output_directory + "/RAxML_info." + file_number)
                        # os.rename("topology_bestTree." + file_number, topology_output_directory + "/Topology_bestTree." + file_number)

                percent_complete += 80 / len(os.listdir(window_directory))
                self.emit(QtCore.SIGNAL('RAX_PER'), percent_complete)

    def run(self):
        try:
            self.window_splitter(self.inputFilename, self.windowSize, self.windowOffset)
        except IOError:
            self.emit(QtCore.SIGNAL('INVALID_ALIGNMENT_FILE'), self.inputFilename)
            return

        self.raxml_windows(self.numBootstraps, self.model, rooted=self.rooted, outgroup=self.outGroup)
        self.emit(QtCore.SIGNAL('RAX_COMPLETE'), None)

        print 'Alignment:', self.inputFilename
        print 'Window Size:', self.windowSize
        print 'Window Offset:', self.windowOffset
        print '# Bootstraps:', self.numBootstraps


if __name__ == '__main__':
    inputFile = "../RAxML_SpeciesTree/RAxML_ST_bestTree.txt"
    windowSize = 500000
    windowOffset = 500000
    numBootstraps = 2

    ro = RAxMLOperations(inputFile, windowSize, windowOffset, numBootstraps=2)

    windows_dir = ro.window_splitter(ro.inputFilename, ro.windowSize, ro.windowOffset)
    ro.raxml_windows()
