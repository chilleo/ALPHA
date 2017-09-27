import subprocess
import os
from natsort import natsorted
import re
import dendropy
from dendropy import Tree
import math
from dendropy.calculate import treecompare
import matplotlib.pyplot as plt
import numpy as np
from PyQt4 import QtCore

"""
Functions:
    __init__(self, output_directory='RAxML_Files', parent=None)
    newick_reformat(self, newick)
    calculate_p_of_gt_given_st(self, species_tree, gene_tree)
    calculate_windows_to_p_gtst(self, species_tree)
    calculate_robinson_foulds(self, species_tree, gene_tree, weighted)
    calculate_windows_to_rf(self, species_tree, weighted)
    stat_scatter(self, stat_map, name, title, xlabel, ylabel)
    calculate_d(self, alignment, window_size, window_offset, taxon1, taxon2, taxon3, taxon4)
    run(self)
~
Chabrielle Allen
Travis Benedict
Peter Dulworth
"""


class StatisticsCalculations(QtCore.QThread):
    def __init__(self, output_directory='RAxML_Files', parent=None):
        super(StatisticsCalculations, self).__init__(parent)
        self.output_directory = output_directory

    def newick_reformat(self, newick):
        """
        Reformat the inputted newick string to work with the PhyloNet jar file
        "(a:2.5,(b:1.0,c:1.0):1.5)" This format works
        "(a:2.0,(b:1.0,c:1.0):1.0);" This format works
        "(a:2.0,(b:1.0,c:1.0)):1.0;" THIS FORMAT DOES NOT WORK --- trees from RAxML are in this format
        Inputs:
        newick --- an incorrectly formatted newick string
        Output:
        newick --- a correctly formatted newick string
        """

        # Find root length and remove it
        pattern = "(?!.*\))(.*?)(?=\;)"

        newick = re.sub(pattern, '', newick)

        return newick

    def calculate_p_of_gt_given_st(self, species_tree, gene_tree):
        """
            Computes the probability that a gene tree occurs given a species tree. If the taxon names between the two trees are not the
            same then the probability returned is 0.0. If trees are the exact same then probability is 1.0.

            Inputs:
                species_tree --- a newick string containing a species tree with branch lengths as outputted by RAxML or inputted by user
                gene_tree --- a newick string containing a gene tree with branch lengths as outputted by RAxML run on windows
            Output:
                p_of_gt_given_st --- the probability that a gene tree occurs given a species tree.
        """

        # If species_tree input is a file read in the newick string
        if os.path.isfile(species_tree):
            with open(species_tree) as f:
                species_tree = f.readline()

        # Check if the species tree is formatted correctly for PhyloNet if not reformat it
        if species_tree[-2] != ")" or species_tree[-1] != ")":
            # species_tree = newick_reformat(species_tree).replace("\n","")
            species_tree = species_tree.replace("\n", "")

        # If gene_tree input is a file read in the newick string
        if os.path.isfile(gene_tree):
            with open(gene_tree) as f:
                gene_tree = f.readline()

        # Check if the gene tree is formatted correctly for PhyloNet if not reformat it
        if gene_tree[-2] != ")" or gene_tree[-1] != ")":
            # gene_tree = newick_reformat(gene_tree).replace("\n","")
            gene_tree = gene_tree.replace("\n", "")

        # IF YOU COMMENT THIS OUT AGAIN EVERYTHING WILL BREAK
        # add quotes to the strings
        species_tree = str(species_tree)
        species_tree = "'" + species_tree + "'"
        gene_tree = "'" + gene_tree + "'"

        # Run PhyloNet jar file
        p = subprocess.Popen("java -jar ./unstable.jar {0} {1}".format(species_tree, gene_tree), stdout=subprocess.PIPE, shell=True)

        # Read output and convert to float
        p_of_gt_given_st = p.stdout.readline()

        return p_of_gt_given_st

    def calculate_windows_to_p_gtst(self, species_tree):
        """
        Calculate p(gt|st) for each window and create a mapping
        of window numbers to probabilities.

        Inputs:
            species_tree --- a newick string containing a species tree with branch lengths as outputted by RAxML or inputted by user
        Output:
            windows_to_p_gtst --- a mapping of window numbers to their p(gt|st).
        """

        # Initialize a mapping
        windows_to_p_gtst = {}

        # Iterate over each folder in the given directory
        for filename in natsorted(os.listdir(self.output_directory)):

            # If file is the file with the best tree newick string
            if os.path.splitext(filename)[0] == "RAxML_bestTree":
                window_num = (os.path.splitext(filename)[1]).replace(".", "")

                gene_tree_filename = os.path.join(self.output_directory, filename)

                p_gtst = self.calculate_p_of_gt_given_st(species_tree, gene_tree_filename)

                # Reformat output
                p_gtst = float(p_gtst.replace('\r', '').replace('\n', ''))

                windows_to_p_gtst[window_num] = p_gtst

        return windows_to_p_gtst

    def calculate_robinson_foulds(self, species_tree, gene_tree, weighted):
        """
        Calculates the Robinson Foulds distances for weighted and unweighted
        trees.

        Input:
        species_tree -- newick file or newick string containing the species tree
        gene_tree   -- newick file or newick string containing the tree to
                          be compared to the species tree
        weighted       -- boolean parameter for whether the files have weights

        Returns:
        The weighted and/or unweighted Robinson Foulds distance of the species
        tree and input tree.
        """

        # taxon names
        tns = dendropy.TaxonNamespace()

        # Create dendropy tree from species tree input file
        if os.path.isfile(species_tree):
            species_tree = Tree.get_from_path(species_tree, 'newick', taxon_namespace=tns)

        # Create dendropy tree from species tree input newick string
        else:
            species_tree = Tree.get_from_string(species_tree, 'newick', taxon_namespace=tns)

        # Create dendropy tree from gene tree input file
        if os.path.isfile(gene_tree):
            gene_tree = Tree.get_from_path(gene_tree, 'newick', taxon_namespace=tns)

        # Create dendropy tree from gene tree input newick string
        else:
            gene_tree = Tree.get_from_string(gene_tree, 'newick', taxon_namespace=tns)

        # both weighted and unweighted foulds distance
        if weighted:
            return treecompare.weighted_robinson_foulds_distance(species_tree, gene_tree), \
                   treecompare.unweighted_robinson_foulds_distance(species_tree, gene_tree)

        # only unweighted foulds distance
        else:
            return treecompare.unweighted_robinson_foulds_distance(species_tree, gene_tree)

    def calculate_windows_to_rf(self, species_tree, weighted):
        """
        Calculate Robinson-Foulds distance for each window and create a
        mapping of window numbers to RF distance.
        Inputs:
        species_tree --- a newick string containing a species tree with
                         branch lengths as outputted by RAxML or inputted
                         by user
        weighted --- a boolean corresponding to calculating the weighted
                     or unweighted RF distance
        Output:
        windows_to_rf --- a mapping of window numbers to their RF distance.
        """

        # Initialize a mapping for the weighted and unweighted RF distance
        windows_to_w_rf = {}
        windows_to_uw_rf = {}

        # Iterate over each folder in the given directory
        for filename in natsorted(os.listdir(self.output_directory)):

            # If file is the file with the best tree newick string
            if os.path.splitext(filename)[0] == "RAxML_bestTree":
                # makes file and calculates rf distance
                window_num = (os.path.splitext(filename)[1]).replace(".", "")

                gene_tree_filename = os.path.join(self.output_directory, filename)

                rf_distance = self.calculate_robinson_foulds(species_tree, gene_tree_filename, weighted)

                if weighted:

                    # Weighted RF
                    windows_to_w_rf[window_num] = rf_distance[0]
                    # Unweighted RF
                    windows_to_uw_rf[window_num] = rf_distance[1]

                else:

                    # Unweighted RF
                    windows_to_uw_rf[window_num] = rf_distance

        # returns weighted and/or unweighted Robinson Foulds mappings
        if weighted:
            return windows_to_w_rf, windows_to_uw_rf

        else:
            return windows_to_uw_rf

    def stat_scatter(self, stat_map, name, title, xlabel, ylabel):
        """
        Creates a scatter plot with the x-axis being the
        windows and the y-axis being the statistic to
        be graphed.

        Input:
        stat_map -- a mapping outputted by either
                    calculate_windows_to_p_gtst or
                    calculate_windows_to_rf ([0] or [1])
        name -- the name of the save file
        title -- the title of the plot
        xlabel -- the label for the x axis
        ylabel -- the label for the y axis

        Returns:
        A scatter plot with windows as the x-axis and
        a statistic as the y-axis.
        """
        # sizes plot circles
        area = math.pi * (3) ** 2

        x_list = []

        # makes x values integers
        xlist = stat_map.keys()
        for j in range(len(xlist)):
            x_list.append(int(xlist[j]))

        x = np.array(x_list)

        # gets y values from dictionary
        ylist = stat_map.values()
        y = np.array(ylist)

        plt.scatter(x, y, s=area, c='#000000', alpha=1)

        # label the axes
        plt.xlabel(xlabel, fontsize=10)
        plt.ylabel(ylabel, fontsize=10)

        plt.title(title, fontsize=15)
        plt.tight_layout()
        plt.savefig(name)

        plt.clf()

    def figureBarPlot(self, data, name, title, labelHeights=False, legend=False, legendNames=(), xTicks=False, groupLabels=()):
        """
            generates bar chart
        """

        figure = plt.figure()

        numberOfBars = len(data)
        ind = np.arange(numberOfBars)  # the x locations for the groups
        width = .667  # the width of the bars
        ax = figure.add_subplot(111)
        colors = [(43.0/255.0, 130.0/255.0, 188.0/255.0), (141.0/255.0, 186.0/255.0, 87.0/255.0), (26.0/255.0, 168.0/255.0, 192.0/255.0), (83.5/255.0, 116.5/255.0, 44.5/255.0)]

        ax.bar(ind, data, width, color=colors)

        plt.title(title, fontsize=15)
        plt.savefig(name)
        plt.show()

        # plt.clf()

    def barPlot(self, data, name, title, xLabel, yLabel, labelHeights=False, legend=False, legendNames=(), xTicks=False, groupLabels=()):
        """
            generates bar chart
        """

        def autoLabel(rects, ax):
            """
            Attach a text label above each bar displaying its height
            """
            for rect in rects:
                height = rect.get_height()
                ax.text(rect.get_x() + rect.get_width() / 2., 1.05 * height, '%d' % int(height), ha='center', va='bottom')


        plt.figure()

        numberOfBars = len(data)
        ind = np.arange(numberOfBars)  # the x locations for the groups
        width = .667  # the width of the bars
        fig, ax = plt.subplots()
        colors = [(43.0/255.0, 130.0/255.0, 188.0/255.0), (141.0/255.0, 186.0/255.0, 87.0/255.0), (26.0/255.0, 168.0/255.0, 192.0/255.0), (83.5/255.0, 116.5/255.0, 44.5/255.0)]
        bars = []

        bars.append(ax.bar(ind, data, width, color=colors))

        ax.set_xticks([])
        if xTicks:
            ax.set_xticks(ind)

        if len(groupLabels) > 0:
            ax.set_xticks(ind)
            ax.set_xticklabels(groupLabels)

        if legend:
            legendBoxes = []
            for bar in bars:
                legendBoxes.append(bar[0])
            ax.legend(legendBoxes, legendNames)

        if labelHeights:
            autoLabel(bars[0], ax)
            autoLabel(bars[1], ax)

        # label the axes
        plt.xlabel(xLabel, fontsize=10)
        plt.ylabel(yLabel, fontsize=10)

        plt.title(title, fontsize=15)
        # plt.tight_layout()
        plt.savefig(name)

        # plt.clf()

    def groupedBarPlot(self, data, name, title, xLabel, yLabel, labelHeights=False, legend=False, legendNames=(), xTicks=False, groupLabels=()):
        """
                    generates bar chart
                """

        def autoLabel(rects, ax):
            """
            Attach a text label above each bar displaying its height
            """
            for rect in rects:
                height = rect.get_height()
                ax.text(rect.get_x() + rect.get_width() / 2., 1.05 * height, '%d' % int(height), ha='center', va='bottom')

        useGroups = True

        print data

        if type(data[0]) == list:
            if len(data[0]) == 1:
                useGroups = False
        else:
            newData = []
            for el in data:
                newData.append([el])
            data = newData
            useGroups = False

        numberOfBarGroups = len(data[0])
        ind = np.arange(numberOfBarGroups)  # the x locations for the groups
        width = (0.6667) / float(len(data))  # the width of the bars
        offset = 0.0
        fig, ax = plt.subplots()
        colors = [(43.0 / 255.0, 130.0 / 255.0, 188.0 / 255.0), (141.0 / 255.0, 186.0 / 255.0, 87.0 / 255.0), (26.0 / 255.0, 168.0 / 255.0, 192.0 / 255.0), (83.5 / 255.0, 116.5 / 255.0, 44.5 / 255.0)]
        groups = []
        padding = 0

        if not useGroups:
            padding = 0.25

        for i in range(len(data)):
            groups.append(ax.bar(ind + offset, data[i], width, color=colors[i % 4]))
            offset += width + padding

        ax.set_xticks([])
        if xTicks:
            ax.set_xticks(ind + ((numberOfBarGroups - 1.0)) * width / 2.0)

        if len(groupLabels) > 0:
            ax.set_xticklabels(groupLabels)

        if legend:
            legendBoxes = []
            for group in groups:
                legendBoxes.append(group[0])
            ax.legend(legendBoxes, legendNames)

        if labelHeights:
            autoLabel(groups[0], ax)
            autoLabel(groups[1], ax)

        # label the axes
        plt.xlabel(xLabel, fontsize=10)
        plt.ylabel(yLabel, fontsize=10)

        plt.title(title, fontsize=15)
        plt.tight_layout()
        plt.savefig(name)
        # plt.show()

        plt.clf()

    def calculate_d(self, alignment, window_size, window_offset, taxon1, taxon2, taxon3, taxon4):
        """
            Calculates the D statistic for the given alignment

            Input:
                i. alignment --- a sequence alignment in phylip format
                ii. window_size --- the size of the desired windows
                iii. window_offset --- the offset that was used to create the windows
                iv. taxon1 --- first taxon for ABBA-BABA test
                v. taxon2 --- second taxon for ABBA-BABA test
                vi. taxon3 --- third taxon for ABBA-BABA test
                vii. taxon4 --- outgroup for ABBA-BABA test

            Output:
                i. d_stat --- the D statistic value
                ii. windows_to_d --- a mapping of window indices to D values
        """

        # Initialize the site index to 0
        site_idx = 0

        # Initialize values for the d statistic numerator and denominator
        d_numerator = 0
        d_denominator = 0

        windows_to_d = {}

        sequence_list = []
        taxon_list = []

        with open(alignment) as f:

            # Create a list of each line in the file
            lines = f.readlines()

            # First line contains the number and length of the sequences
            first_line = lines[0].split()
            number_of_sequences = int(first_line[0])
            length_of_sequences = int(first_line[1])

        # Initialize a count for the total number of windows
        num_windows = 0

        i = 0

        if window_size > length_of_sequences:
            num_windows = 1
            window_size = length_of_sequences
        else:
            # Determine the total number of windows needed
            while (i + window_size - 1 < length_of_sequences):
                i += window_size
                num_windows += 1

        for line in lines[1:]:
            # Add each sequence to a list
            sequence = line.split()[1]
            sequence_list.append(sequence)

            # Add each taxon to a list
            taxon = line.split()[0]
            taxon_list.append(taxon)

        # Get the index of each taxon in the inputted order
        taxon1_idx = taxon_list.index(taxon1)
        taxon2_idx = taxon_list.index(taxon2)
        taxon3_idx = taxon_list.index(taxon3)
        taxon4_idx = taxon_list.index(taxon4)

        # Initialize values for the d statistic numerator and denominator for each window
        numerator_window = 0
        denominator_window = 0

        # Iterate over each window
        for window in range(num_windows):

            # Iterate over the indices in each window
            for window_idx in range(window_size):

                site = []

                # Iterate over each sequence in the alignment
                for sequence in sequence_list:
                    # Add each base in a site to a list
                    site.append(sequence[site_idx])

                # Get the genetic sites in the correct order
                P1, P2, P3, O = site[taxon1_idx], site[taxon2_idx], site[taxon3_idx], site[taxon4_idx]

                # Case of ABBA
                if P1 == O and P2 == P3 and P1 != P2:
                    ABBA = 1
                    BABA = 0

                # Case of BABA
                elif P1 == P3 and P2 == O and P1 != P2:
                    ABBA = 0
                    BABA = 1

                # Neither case
                else:
                    ABBA = 0
                    BABA = 0

                numerator_window += (ABBA - BABA)
                denominator_window += (ABBA + BABA)

                # Increment the site index
                site_idx += 1

                # Handle case for when the final window is shorter than the others
                if site_idx > length_of_sequences:
                    break

            # Calculate d statistic for the window
            d_window = numerator_window / float(denominator_window)

            # Map the window index to its D statistic
            windows_to_d[window] = d_window

            # Add the numerator and denominator of each window to the overall numerator and denominator
            d_numerator += numerator_window
            d_denominator += denominator_window

            # Reset numerator and denominator
            numerator_window = 0
            denominator_window = 0

            # Account for overlapping windows
            site_idx += (window_offset - window_size)

            percent_complete = (float(window + 1) / float(num_windows)) * 100
            self.emit(QtCore.SIGNAL('D_PER'), percent_complete)

        d_stat = d_numerator / float(d_denominator)

        self.emit(QtCore.SIGNAL('D_FINISHED'), d_stat, windows_to_d)

    def run(self):
        try:
            self.calculate_d(self.dAlignment, self.dWindowSize, self.dWindowOffset, self.taxons[0], self.taxons[1], self.taxons[2], self.taxons[3])
        except IOError:
            self.emit(QtCore.SIGNAL('INVALID_ALIGNMENT_FILE'), 'Invalid File.', 'Invalid alignment file. Please choose another.', self.dAlignment)
        except ZeroDivisionError:
            self.emit(QtCore.SIGNAL('INVALID_ALIGNMENT_FILE'), 'Invalid Taxon Selection.', 'Please make sure that all taxons are unique.', self.dAlignment)
        finally:
            return


if __name__ == '__main__':

    sc = StatisticsCalculations()

    sc.output_directory = "../RAxML_Files"

    sc.calculate_windows_to_p_gtst("(C:0.00773900199203547429,(G:0.00922097500041624447,O:0.04766300468995082057):0.00139495391245007404,H:0.00972794777559046753):0.0;")
