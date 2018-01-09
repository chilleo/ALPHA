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
from sys import platform
import itertools
import ete3
from ete3 import Tree
import copy
from collections import defaultdict

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


        if platform == 'darwin':
            # IF YOU COMMENT THIS OUT AGAIN EVERYTHING WILL BREAK
            # add quotes to the strings
            species_tree = str(species_tree)
            species_tree = "'" + species_tree + "'"
            gene_tree = "'" + gene_tree + "'"

        # Run PhyloNet jar file
        dir_path = os.path.dirname(os.path.realpath(__file__))
        j = os.path.join(dir_path, "Unstable.jar")
        p = subprocess.Popen("java -jar {0} {1} {2}".format(j,species_tree, gene_tree), stdout=subprocess.PIPE, shell=True)

        # Read output and convert to float
        (p_of_gt_given_st, err) = p.communicate()

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

            # Handle division by zero appropriately
            if denominator_window != 0:
                # Calculate d statistic for the window
                d_window = numerator_window / float(denominator_window)
            else:
                d_window = 0

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

        if d_denominator != 0:
            d_stat = d_numerator / float(d_denominator)
        else:
            d_stat = 0

        self.emit(QtCore.SIGNAL('D_FINISHED'), d_stat, windows_to_d)

    # Generalized D statistic functions

    def generate_network_tree(self, inheritance, species_tree, network_map):
        """
        Creates a network tree based on the species tree
        and the two leaves to be connected.
        Inputs:
        inheritance  -- inputted tuple containing inheritance
                        probability ex. (0.7, 0.3)
        species_tree -- generated or inputted file or newick
                        string
        network_map  -- inputted mapping of leaves where nodes
                        will be added
        Returns:
        A newick string network with the added nodes.
        """
        # check for a species tree file
        if os.path.isfile(species_tree):
            with open(species_tree) as f:
                s_tree = f.readline()

        # check for a species tree string
        else:
            s_tree = species_tree

        # get taxa for the edge in the network
        start = network_map.keys()[0]
        end = network_map[start]

        # add nodes into tree in proper format
        network = s_tree.replace(start, '((' + start + ')#H1:0::' + str(inheritance[0]) + ')')
        network = network.replace(end, '(#H1:0::' + str(inheritance[1]) + ',' + end + ')')

        return network

    # Generate all unique trees functions

    def gendistinct(self, n):
        """
        Generate all full binary trees with n leaves
        Input:
        n --- the number of leaves
        Output:
        dp[-1] --- the set of all full binary trees with n nodes
        """

        leafnode = '(.)'
        dp = []
        newset = set()
        newset.add(leafnode)
        dp.append(newset)

        for i in range(1, n):
            newset = set()
            for j in range(i):
                for leftchild in dp[j]:
                    for rightchild in dp[i - j - 1]:
                        newset.add('(' + '.' + leftchild + rightchild + ')')
            dp.append(newset)

        return dp[-1]

    def generate_all_trees(self, taxa):
        """
        Create all trees given a set of taxa
        Inputs:
        taxa --- a set of the taxa to be used for leaf names
        Output:
        trees --- the set of all trees over the taxa
        """

        # Regex pattern for identifying leaves next to a clade in newick string
        pattern = "([\)][a-zA-Z0-9_.-])"

        # Generate all distinct binary trees
        trees = self.gendistinct(len(taxa))

        # Get all possible permutations of the taxa
        taxa_orders = itertools.permutations(taxa)
        taxa_orders = list(taxa_orders)

        all_trees = []

        # Iterate over each tree in the set
        for tree in trees:
            # print 'tree', tree
            # Reformat the tree
            tree = tree.replace('.', '')

            # Iterate over each permutation of taxa
            for taxa_perm in taxa_orders:

                # Create a copy of the tree
                bi_tree = tree

                # replace the leaves with taxons and reformat string
                for i in range(len(taxa_perm)):
                    taxon = taxa_perm[i] + ","
                    bi_tree = bi_tree.replace("()", taxon, 1)

                bi_tree = bi_tree.replace(",)", ")")

                # Find all instances of a ")" followed by a taxon and add a "," between
                clades = re.findall(pattern, bi_tree)
                for clade in clades:
                    taxon = clade[1]
                    bi_tree = bi_tree.replace(clade, ")," + taxon)
                bi_tree = bi_tree.replace(")(", "),(")
                bi_tree = bi_tree + ";"
                all_trees.append(bi_tree)

                # print bi_tree

        return all_trees

    def generate_unique_trees(self, taxa, outgroup):
        """
        Generate the set of unique trees over a set of taxa with an outgroup
        Inputs:
        taxa --- a list of taxa to be used as the leaves of trees
        outgroup --- the outgroup to root at
        Output:
        unique_newicks --- a set of all unique topologies over the given taxa
        """

        # Regular expression for identifying floats
        float_pattern = "([+-]?\\d*\\.\\d+)(?![-+0-9\\.])"
        # Regular expressions for removing branch lengths and confidence values
        pattern2 = "([\:][\\d])"
        pattern3 = "([\)][\\d])"

        # Create a set for unique trees
        unique_trees = set([])

        unique_newicks = set([])

        all_trees = self.generate_all_trees(taxa)

        # Iterate over each tree in all_trees
        for tree in all_trees:

            tree = Tree(tree)
            tree.set_outgroup(outgroup)

            is_unique = True

            # Iterate the unique trees for comparison
            for unique_tree in unique_trees:

                # Compute robinson-foulds distance
                rf_distance = tree.robinson_foulds(unique_tree)[0]

                # If rf distance is 0 the tree is not unique
                if rf_distance == 0:
                    is_unique = False

            if is_unique:
                unique_trees.add(tree)

        # Iterate over the trees
        for tree in unique_trees:
            # Get newick strings from the tree objects
            tree = tree.write()

            # Get rid of branch lengths in the newick strings
            tree = (re.sub(float_pattern, '', tree))
            tree = (re.sub(pattern2, '', tree)).replace(":", "")
            tree = (re.sub(pattern3, ')', tree))

            tree = self.outgroup_reformat(tree, outgroup)

            # Add the newick strings to the set of unique newick strings
            unique_newicks.add(tree)

        return unique_newicks

    # L-Statistic Calculations Functions

    def outgroup_removal(self, newick, outgroup):
        """
        Move the location of the outgroup in a newick string to be at the end of the string
        Inputs:
        newick --- a newick string to be reformatted
        outgroup --- the outgroup
        """

        # Replace the outgroup and comma with an empty string
        newick = newick.replace("," + outgroup, "")

        newick = newick[1:-2] + ";"

        return newick

    def calculate_newicks_to_stats(self, species_tree, species_network, unique_trees, outgroup):
        """
        Compute p(g|S) and p(g|N) for each g in unique_trees and 
        map the tree newick string to those values
        Inputs:
        species_tree --- the species tree newick string for the taxa
        species_network --- the network newick string derived from adding a branch to the species tree between the interested taxa
        unique_trees --- the set of all unique topologies over n taxa
        outgroup --- the outgroup
        Output:
        trees_to_pgS--- a mapping of tree newick strings to their p(g|S) values 
        trees_to_pgN--- a mapping of tree newick strings to their p(g|N) values
        """

        trees_to_pgS = {}
        trees_to_pgN = {}
        trees_to_pgS_noO = {}
        trees_to_pgN_noO = {}

        species_tree_noO = self.outgroup_removal(species_tree, outgroup)
        species_network_noO = self.outgroup_removal(species_network, outgroup)

        # Iterate over the trees
        for tree in unique_trees:
            # Run PhyloNet p(g|S) jar file
            p = subprocess.Popen("java -jar ../unstable.jar {0} {1}".format(species_tree, tree), stdout=subprocess.PIPE,
                                 shell=True)

            # Read output and convert to float
            p_of_g_given_s = float(p.stdout.readline())

            # Run PhyloNet p(g|N) jar file
            p = subprocess.Popen("java -jar ../unstable.jar {0} {1}".format(species_network, tree),
                                 stdout=subprocess.PIPE,
                                 shell=True)

            # Read output and convert to float
            p_of_g_given_n = float(p.stdout.readline())

            # Calculate for trees without outgroup
            tree_noO = self.outgroup_removal(tree, outgroup)

            # Run PhyloNet p(g|S) jar file
            p = subprocess.Popen("java -jar ../unstable.jar {0} {1}".format(species_tree_noO, tree_noO),
                                 stdout=subprocess.PIPE,
                                 shell=True)

            # Read output and convert to float
            p_of_g_given_s_noO = float(p.stdout.readline())

            # Run PhyloNet p(g|N) jar file
            p = subprocess.Popen("java -jar ../unstable.jar {0} {1}".format(species_network_noO, tree_noO),
                                 stdout=subprocess.PIPE,
                                 shell=True)

            # Read output and convert to float
            p_of_g_given_n_noO = float(p.stdout.readline())

            trees_to_pgS[tree] = p_of_g_given_s
            trees_to_pgN[tree] = p_of_g_given_n
            trees_to_pgS_noO[tree] = p_of_g_given_s_noO
            trees_to_pgN_noO[tree] = p_of_g_given_n_noO

        return trees_to_pgS, trees_to_pgN, trees_to_pgS_noO, trees_to_pgN_noO

    def determine_interesting_trees(self, trees_to_pgS, trees_to_pgN):
        """
        Get the subset of trees who are initially equal based on p(g|S) but unequal based on p(g|N)
        Input:
        trees_to_pgS--- a mapping of tree newick strings to their p(g|S) values 
        trees_to_pgN--- a mapping of tree newick strings to their p(g|N) values
        Output:
        interesting_trees --- the subset of tree topologies to look at for determining introgression
        """

        # Initialize a set to contain all tree that are equal based on p(g|S)
        possible_trees = []

        # Compare the probability of each tree to the probability of every other tree
        for tree1 in trees_to_pgS:

            equal_trees = set([])

            for tree2 in trees_to_pgS:

                if trees_to_pgS[tree1] == trees_to_pgS[tree2]:
                    equal_trees.add(tree2)

            if len(equal_trees) > 1:
                # Add the equal trees to the set of possible trees
                possible_trees.append(equal_trees)

        valuable_trees = []

        # Iterate over each set of equal trees
        for equal_trees in possible_trees:

            unequal_trees = set([])

            # Compare the p(g|N) values
            for tree1 in equal_trees:

                for tree2 in equal_trees:

                    if trees_to_pgN[tree1] != trees_to_pgN[tree2]:
                        unequal_trees.add(tree1)
                        unequal_trees.add(tree2)

            if len(unequal_trees) > 0:
                valuable_trees.append(unequal_trees)

        minimal_size = float("inf")

        # Get the minimal subset of interesting trees
        for trees in valuable_trees:
            if len(trees) < minimal_size:
                minimal_size = len(trees)
                interesting_trees = trees

        return interesting_trees

    # Site Pattern Functions

    def outgroup_reformat(self, newick, outgroup):
        """
        Move the location of the outgroup in a newick string to be at the end of the string
        Inputs:
        newick --- a newick string to be reformatted
        outgroup --- the outgroup
        """

        # Replace the outgroup and comma with an empty string
        newick = newick.replace(outgroup + ",", "")

        newick = newick[:-2] + "," + outgroup + ");"

        return newick

    def pattern_inverter(self, patterns):
        """
        Switches "A"s to "B"s and "B"s to "A" in a site pattern excluding the outgroup
        Inputs:
        patterns --- a list of site patterns
        Output:
        inverted --- a list of the inverted patterns
        """

        inverted = []

        # Iterate over the patterns
        for pattern in patterns:

            a_count = 0
            b_count = 0

            inverted_pattern = []

            # Iterate over each site in the pattern
            for site in pattern:

                if site == "A":
                    inverted_pattern.append("B")
                    b_count += 1

                elif site == "B":
                    inverted_pattern.append("A")
                    a_count += 1

            if inverted_pattern[-1] != "A":
                # Change the last site to an "A"
                inverted_pattern[-1] = "A"
                b_count -= 1
                a_count += 1

            if a_count > 1 and b_count > 0:
                inverted.append(inverted_pattern)

        return inverted

    def pattern_string_generator(self, patterns):
        """
        Creates a list of viable pattern strings that are easier to read
        Input:
        patterns --- a list of lists of individual characters e.g. [["A","B","B","A"],["B","A","B","A"]]
        Output:
        pattern_strings --- a list of lists of strings e.g. [["ABBA"],["BABA"]]
        """

        # Convert the site pattern lists to strings
        pattern_strings = []
        while patterns:

            a_count = 0
            b_count = 0
            pattern_str = ""
            pattern = patterns.pop()

            for site in pattern:

                if site == "A":
                    b_count += 1

                elif site == "B":
                    a_count += 1

                pattern_str += site

            if a_count > 0 and b_count > 0:
                pattern_strings.append(pattern_str)

        return pattern_strings

    def site_pattern_generator(self, taxa_order, newick, outgroup):
        """
        Generate the appropriate AB list patterns
        Inputs:
        taxa_order --- the desired order of the taxa
        newick --- the newick string to generate site patterns for
        outgroup --- the outgroup of the tree
        Output:
        finished_patterns --- the list of site patterns generated for the newick string
        """
        # Create a tree object
        tree = ete3.Tree(newick, format=1)

        # Initialize containers for the final patterns and patterns being altered
        final_site_patterns = []

        # Keep a count of clades in the tree that contain 2 leaves
        clade_count = 0

        # Number of possible patterns is number of taxa - 2 + the number of clades
        num_patterns = len(taxa_order) - 2

        # Initialize pattern to be a list of strings
        pattern = ["B" for x in range(len(taxa_order))]

        # Create list of nodes in order of appearance
        nodes = []
        for node in tree.traverse("postorder"):
            # Add node name to list of nodes
            nodes.append(node.name)

        nodes = list(reversed(nodes))

        if nodes[2] == "" and nodes[3] == "":
            nodes = []
            for node in tree.traverse("preorder"):
                # Add node name to list of nodes
                nodes.append(node.name)

            nodes = list(reversed(nodes))

        # Keep track of visited leaves
        seen_leaves = []

        # Create a list of patterns that are duplicates
        duplicates = []

        # Iterate over the order that the nodes occur beginning at the root
        for node_idx in range(len(nodes)):

            node = nodes[node_idx]

            # If the node is the outgroup add A to the end of the pattern
            if node == outgroup:
                pattern[-1] = "A"
                # Add outgroup to the seen leaves
                seen_leaves.append(node)

            # Else if the node is a leaf and is adjacent to the outgroup
            elif node != "" and seen_leaves[-1] == outgroup and outgroup in seen_leaves:

                # If the next node is a leaf a clade has been found
                if nodes[node_idx + 1] != "":
                    node2 = nodes[node_idx + 1]

                    # Get the indices of the leaves in the pattern
                    pat_idx1 = taxa_order.index(node)
                    pat_idx2 = taxa_order.index(node2)

                    # Set those pattern indices to "A"
                    pattern[pat_idx1] = "A"
                    pattern[pat_idx2] = "A"

                    clade_count += 1

                    # If there is a clade besides the first one then duplicate it in the list
                    final_site_patterns.append(pattern)
                    duplicates.append(pattern)

                    seen_leaves.append(node)
                    seen_leaves.append(node2)

                    # Get the index that final clade occurs at
                    end_idx = node_idx + 1
                    break

                # Otherwise there is no clade
                else:
                    # Get the index of the leaf in the pattern
                    pat_idx = taxa.index(node)

                    # Set those pattern indices to "A"
                    pattern[pat_idx] = "A"

                    seen_leaves.append(node)

                    # Get the index that final leaf occurs at
                    end_idx = node_idx
                    break

        num_patterns = num_patterns + clade_count
        # All patterns can be derived from the pattern with the most B's
        working_patterns = [pattern for x in range(num_patterns)]

        # Pop a pattern off of working patterns and add it to the final site patterns
        final_site_patterns.append(working_patterns.pop())

        # Iterate over each pattern in working patterns and change them
        while working_patterns:

            # Get a pattern and copy it
            pattern = copy.deepcopy(working_patterns.pop())

            # Iterate over the order that the nodes occur beginning at the last clade or leaf
            for node_idx in range(end_idx + 1, len(nodes)):

                # If the last clade is reached break
                if node_idx == len(nodes) - 1:

                    if node != "":
                        # Get the index of the leaf in the pattern
                        pat_idx1 = taxa_order.index(node)

                        # Set those pattern indices to "A"
                        pattern[pat_idx1] = "A"

                        # Get the index that final leaf occurs at
                        end_idx = node_idx
                        break

                    else:
                        break

                node = nodes[node_idx]

                # If the next node is a leaf a clade has been found
                if node != "" and nodes[node_idx + 1] != "":
                    node2 = nodes[node_idx + 1]

                    # Get the indices of the leaves in the pattern
                    pat_idx1 = taxa_order.index(node)
                    pat_idx2 = taxa_order.index(node2)

                    # Set those pattern indices to "A"
                    pattern[pat_idx1] = "A"
                    pattern[pat_idx2] = "A"

                    clade_count += 1

                    # If there is a clade besides the first one then duplicate it in the list
                    final_site_patterns.append(pattern)
                    duplicates.append(pattern)

                    # Get the index that final clade occurs at
                    end_idx = node_idx + 1
                    break

                # Else if the node is a leaf
                elif node != "":
                    # Get the index of the leaf in the pattern
                    pat_idx1 = taxa_order.index(node)

                    # Set those pattern indices to "A"
                    pattern[pat_idx1] = "A"

                    # Get the index that final leaf occurs at
                    end_idx = node_idx
                    break

            # Add the altered pattern to the final site patterns
            final_site_patterns.append(pattern)
            duplicates.append(pattern)

            # Update the working patterns to be the same as the most recent pattern
            working_patterns = [pattern for x in range(num_patterns - len(final_site_patterns))]

        # Create a list of patterns without duplicates
        finished_patterns = []

        # Iterate over each pattern and determine which ones are duplicates
        for pattern in final_site_patterns:

            if pattern not in finished_patterns:
                finished_patterns.append(pattern)

            else:
                duplicates.append(pattern)

        # This may need to change double check with Chill Leo on this
        # if clade_count > 1:

        duplicates = finished_patterns

        # Invert all duplicate patterns
        inverted_patterns = self.pattern_inverter(duplicates)

        # Iterate over the inverted patterns and add them to finished patterns
        for pattern in inverted_patterns:

            if pattern not in finished_patterns:
                finished_patterns.append(pattern)

        finished_patterns = self.pattern_string_generator(finished_patterns)

        return finished_patterns

    def newicks_to_patterns_generator(self, taxa_order, newicks):
        """
        Generate the site patterns for each newick string and map the strings to their patterns
        Inputs:
        taxa_order --- the desired order of the taxa
        newicks --- a list of newick strings
        Output:
        newicks_to_patterns --- a mapping of newick strings to their site patterns
        """

        # Determine the outgroup of the tree
        outgroup = taxa_order[-1]

        newicks_to_patterns = {}

        # Iterate over the newick strings
        for newick in newicks:
            newicks_to_patterns[newick] = self.site_pattern_generator(taxa_order, newick, outgroup)

        return newicks_to_patterns

    # Interesting sites functions

    def calculate_pattern_probabilities(self, newicks_to_patterns, newicks_to_pgS, newicks_to_pgN):
        """
        Creates a mapping of site patterns to their total p(g|S) values across all gene trees and 
        a mapping of site patterns to their total p(g|N) values across all gene trees
        Inputs:
        newicks_to_patterns --- a mapping of tree newick strings to their site patterns
        newicks_to_pgS--- a mapping of tree newick strings to their p(g|S) values 
        newicks_to_pgN--- a mapping of tree newick strings to their p(g|N) values
        Outputs:
        patterns_to_pgS --- a mapping of site patterns to their total p(g|S) value
        patterns_to_pgN --- a mapping of site patterns to their total p(g|N) value
        """

        patterns_to_pgS = {}
        patterns_to_pgN = {}

        # Iterate over each newick string
        for newick in newicks_to_patterns:
            # Iterate over each site pattern of a tree
            for pattern in newicks_to_patterns[newick]:

                # Initialize a probability for each pattern if it does not have one
                if pattern not in patterns_to_pgS:
                    patterns_to_pgS[pattern] = newicks_to_pgS[newick]
                    patterns_to_pgN[pattern] = newicks_to_pgN[newick]

                # Otherwise add to the existing probability
                else:
                    patterns_to_pgS[pattern] += newicks_to_pgS[newick]
                    patterns_to_pgN[pattern] += newicks_to_pgN[newick]

        return patterns_to_pgS, patterns_to_pgN

    def determine_patterns(self, patterns_to_pgS, patterns_to_pgN1, patterns_to_pgN2):
        """
        Determine which patterns are useful in determining introgression
        Inputs:
        patterns_to_pgS --- a mapping of site patterns to their total p(g|S) value
        patterns_to_pgN1 --- a mapping of site patterns to their total p(g|N) value for a network
        patterns_to_pgN2 --- a mapping of site patterns to their total p(g|N) value for a different network
        Outputs:
        terms1 --- a set of patterns to count and add to each other to determine introgression
        terms2 --- a set of other patterns to count and add to each other to determine introgression
        """

        # Initialize sets for the patterns of interest
        interesting_patterns = set([])
        terms1 = set([])
        terms2 = set([])

        # Iterate over each pattern to determine the patterns of interest
        for pattern in patterns_to_pgS:

            tree_probability = patterns_to_pgS[pattern]

            # If either network probability is not equal to the tree probability and the network probabilities are not equal the pattern is of interest
            if (patterns_to_pgN1[pattern] != tree_probability or patterns_to_pgN2[pattern] != tree_probability) and \
                            patterns_to_pgN1[pattern] != patterns_to_pgN2[pattern]:
                interesting_patterns.add(pattern)

        # Iterate over each interesting pattern and determine which set of terms to add the pattern to
        for pattern in interesting_patterns:

            if patterns_to_pgN1[pattern] > patterns_to_pgN2[pattern]:
                terms1.add(pattern)

            elif patterns_to_pgN2[pattern] > patterns_to_pgN1[pattern]:
                terms2.add(pattern)

        return terms1, terms2

    def generate_statistic_string(self, patterns_of_interest):
        """
        Create a string representing the statistic for determining introgression like "(ABBA - BABA)/(ABBA + BABA)"
        Input:
        patterns_of_interest --- a tuple containing the sets of patterns used for determining a statistic
        Output:
        L_statistic --- a string representation of the statistic
        """

        calculation = []

        # Iterate over each set of patterns
        for pattern_set in patterns_of_interest:
            term = "("

            # Combine each term with a "+"
            for pattern in pattern_set:
                term = term + pattern + " + "
            term = term[:-3] + ")"
            calculation.append(term)

        L_statistic = "({0} - {1}) / ({0} + {1})".format(calculation[0], calculation[1])

        return L_statistic

    # Function for calculating statistic

    def calculate_L(self, alignment, taxa_order, patterns_of_interest):
        """
        Calculates the L statistic for the given alignment
        Input:
        alignment --- a sequence alignment in phylip format
        taxa_order --- the desired order of the taxa
        patterns_of_interest --- a tuple containing the sets of patterns used for determining a statistic
        Output:
        l_stat --- the L statistic value
        """

        # Separate the patterns of interest into their two terms
        terms1 = patterns_of_interest[0]
        terms2 = patterns_of_interest[1]

        terms1_counts = defaultdict(int)
        terms2_counts = defaultdict(int)

        sequence_list = []
        taxon_list = []

        with open(alignment) as f:

            # Create a list of each line in the file
            lines = f.readlines()

            # First line contains the number and length of the sequences
            first_line = lines[0].split()
            length_of_sequences = int(first_line[1])

        for line in lines[1:]:
            # Add each sequence to a list
            sequence = line.split()[1]
            sequence_list.append(sequence)

            # Add each taxon to a list
            taxon = line.split()[0]
            taxon_list.append(taxon)

        # The outgroup is the last taxa in taxa order
        outgroup = taxa_order[-1]

        # Iterate over the site indices
        for site_idx in range(length_of_sequences):

            # Map each taxa to the base at a given site
            taxa_to_site = {}

            # Create a set of the bases at a given site to determine if the site is biallelic
            bases = set([])

            # Iterate over each sequence in the alignment
            for sequence, taxon in zip(sequence_list, taxon_list):
                # Map each taxon to the corresponding base at the site
                base = sequence[site_idx]
                taxa_to_site[taxon] = base
                bases.add(base)

            if len(bases) == 2:

                # Create the pattern that each site has
                site_pattern = []

                # The ancestral gene is always the same as the outgroup
                ancestral = taxa_to_site[outgroup]

                # Iterate over each taxon
                for taxon in taxa_order:
                    nucleotide = taxa_to_site[taxon]

                    # Determine if the correct derived/ancestral status of each nucleotide
                    if nucleotide == ancestral:
                        site_pattern.append("A")
                    else:
                        site_pattern.append("B")

                # Convert the site pattern to a string
                site_string = self.pattern_string_generator([site_pattern])[0]

                # If the site string is a pattern of interest add to its count for one of the terms
                if site_string in terms1:
                    terms1_counts[site_string] += 1

                elif site_string in terms2:
                    terms2_counts[site_string] += 1

        terms1_total = sum(terms1_counts.values())
        terms2_total = sum(terms2_counts.values())

        numerator = terms1_total - terms2_total
        denominator = terms1_total + terms2_total

        if denominator != 0:
            l_stat = numerator / float(denominator)
        else:
            l_stat = 0

        return l_stat

    def l_statistic(self, alignment, taxa, species_tree, reticulations):
        """
        Calculates the L statistic for the given alignment
        Input:
        alignment --- a sequence alignment in phylip format
        taxa --- a list of the taxa in the desired order
        species_tree --- the inputted species tree over the given taxa
        reticulations a tuple containing two dictionaries mapping the start leaves to end leaves
        Output:
        l_stat --- the L statistic value
        """

        # The outgroup is the last taxon in the list of taxa
        outgroup = taxa[-1]

        # Generate all unique trees over the given topology
        unique = self.generate_unique_trees(taxa, outgroup)

        # Map the tree newick strings to their site patterns
        newick_patterns = self.newicks_to_patterns_generator(taxa, unique)

        # Create species networks
        network_map1, network_map2 = reticulations[0], reticulations[1]
        network1 = self.generate_network_tree((0.3, 0.7), species_tree, network_map1)
        network2 = self.generate_network_tree((0.3, 0.7), species_tree, network_map2)

        trees_to_pgS, trees_to_pgN, trees_to_pgS_noO, trees_to_pgN_noO = self.calculate_newicks_to_stats(species_tree,
                                                                                                    network1,
                                                                                                    unique, outgroup)
        patterns_pgS, patterns_pgN1 = self.calculate_pattern_probabilities(newick_patterns, trees_to_pgS, trees_to_pgN)

        trees_to_pgS, trees_to_pgN, trees_to_pgS_noO, trees_to_pgN_noO = self.calculate_newicks_to_stats(species_tree,
                                                                                                    network2,
                                                                                                    unique, outgroup)
        patterns_pgS, patterns_pgN2 = self.calculate_pattern_probabilities(newick_patterns, trees_to_pgS, trees_to_pgN)

        patterns_of_interest = self.determine_patterns(patterns_pgS, patterns_pgN1, patterns_pgN2)

        l_stat = self.calculate_L(alignment, taxa, patterns_of_interest)

        # return l_stat
        self.emit(QtCore.SIGNAL('L_Finished'), l_stat)

    # alignment = "C:\\Users\\travi\\Documents\\PhyloVis\\testFiles\\ChillLeo-Copy.phylip"
    # taxa = ["P1", "P2", "P3", "O"]
    # species_tree = "(((P1:0.8,P2:0.8):0.8,P3:0.8),O);"
    # reticulations = ({"P3": "P2"}, {"P3": "P1"})
    #
    # print l_statistic(alignment, taxa, species_tree, reticulations)

    def run(self):
        try:
            self.calculate_d(self.dAlignment, self.dWindowSize,
                             self.dWindowOffset, self.taxons[0], self.taxons[1], self.taxons[2], self.taxons[3])
        except IOError:
            self.emit(QtCore.SIGNAL('INVALID_ALIGNMENT_FILE'),
                      'Invalid File.', 'Invalid alignment file. Please choose another.', self.dAlignment)
        except ZeroDivisionError:
            self.emit(QtCore.SIGNAL('INVALID_ALIGNMENT_FILE'),
                      'Invalid Taxon Selection.', 'Please make sure that all taxons are unique.', self.dAlignment)
        finally:
            return


if __name__ == '__main__':

    sc = StatisticsCalculations()

    sc.output_directory = "../RAxML_Files"

    sc.calculate_windows_to_p_gtst("(C:0.00773900199203547429,(G:0.00922097500041624447,O:0.04766300468995082057):0.00139495391245007404,H:0.00972794777559046753):0.0;")
