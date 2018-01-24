import re
import os
import itertools
import ete3
import copy
import subprocess
from collections import defaultdict
from sys import platform
from scipy import stats
from ete3 import Tree
from natsort import natsorted
from Bio import AlignIO

"""
Functions:

~
Chabrielle Allen
Travis Benedict
Peter Dulworth
"""


def generate_network_tree(inheritance, species_tree, reticulations):
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
            network = f.readline()

    # check for a species tree string
    else:
        network = species_tree

    for i in range(len(reticulations)):
        # get taxa for the edge in the network
        start = reticulations[i][0]
        end = reticulations[i][1]

        # add nodes into tree in proper format
        network = network.replace(start, '((' + start + ')#H' + str(i+1) + ':0::' + str(inheritance[0]) + ')')
        network = network.replace(end, '(#H' + str(i+1) + ':0::' + str(inheritance[1]) + ',' + end + ')')

    return network


##### Generate all unique trees functions


def genDistinct(n):
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


def generate_all_trees(taxa):
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
    trees = genDistinct(len(taxa))

    # Get all possible permutations of the taxa
    taxa_orders = itertools.permutations(taxa)
    taxa_orders = list(taxa_orders)

    all_trees = []

    # Iterate over each tree in the set
    for tree in trees:
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

    return all_trees


def generate_unique_trees(taxa, outgroup):
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

    all_trees = generate_all_trees(taxa)

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

        tree = outgroup_reformat(tree, outgroup)

        # Add the newick strings to the set of unique newick strings
        unique_newicks.add(tree)

    return unique_newicks


###### Statistics Calculations Functions


def calculate_pgtst(species_tree, gene_tree):
    """
    Calculate p(gt|st) or p(gt|sn)
    Input:
    species_tree --- a species tree or network in newick format
    gene_tree --- a gene tree in newick format
    Output:
    pgtst --- p(gt|st) or p(gt|sn)
    """

    # Run PhyloNet p(g|S) jar file
    p = subprocess.Popen("java -jar unstable.jar {0} {1}".format(species_tree, gene_tree), stdout=subprocess.PIPE,
                         shell=True)

    # Read output and convert to float
    pgtst = float(p.stdout.readline())

    return pgtst


def calculate_newicks_to_stats(species_tree, species_network, unique_trees):
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


    if platform == 'darwin':
        # IF YOU COMMENT THIS OUT AGAIN EVERYTHING WILL BREAK
        # add quotes to the strings
        species_tree = str(species_tree)
        species_tree = "'" + species_tree + "'"
        species_network = str(species_network)
        species_network = "'" + species_network + "'"

    # Iterate over the trees
    for tree in unique_trees:
        if platform == 'darwin':
            # IF YOU COMMENT THIS OUT AGAIN EVERYTHING WILL BREAK
            # add quotes to the strings
            tree = "'" + tree + "'"

        dir_path = os.path.dirname(os.path.realpath(__file__))
        j = os.path.join(dir_path, "Unstable.jar")

        # Run PhyloNet p(g|S) jar file
        p = subprocess.Popen("java -jar {0} {1} {2}".format(j, species_tree, tree), stdout=subprocess.PIPE,
                             shell=True)

        # Read output and convert to float
        p_of_g_given_s = float(p.stdout.readline())

        # Run PhyloNet p(g|N) jar file
        p = subprocess.Popen("java -jar {0} {1} {2}".format(j, species_network, tree), stdout=subprocess.PIPE,
                             shell=True)

        # Read output and convert to float
        p_of_g_given_n = float(p.stdout.readline())

        trees_to_pgS[tree] = p_of_g_given_s
        trees_to_pgN[tree] = p_of_g_given_n

    return trees_to_pgS, trees_to_pgN#, trees_to_pgS_noO, trees_to_pgN_noO

def determine_interesting_trees(trees_to_pgS, trees_to_pgN):
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


##### Site Pattern Functions


def outgroup_reformat(newick, outgroup):
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


def pattern_inverter(patterns):
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


def pattern_string_generator(patterns):
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


def site_pattern_generator(taxa_order, newick, outgroup):
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
                pat_idx = taxa_order.index(node)

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
    inverted_patterns = pattern_inverter(duplicates)

    # Iterate over the inverted patterns and add them to finished patterns
    for pattern in inverted_patterns:

        if pattern not in finished_patterns:
            finished_patterns.append(pattern)

    finished_patterns = pattern_string_generator(finished_patterns)

    return finished_patterns


def newicks_to_patterns_generator(taxa_order, newicks):
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
        newicks_to_patterns[newick] = site_pattern_generator(taxa_order, newick, outgroup)

    return newicks_to_patterns


##### Interesting sites functions


def calculate_pattern_probabilities(newicks_to_patterns, newicks_to_pgS, newicks_to_pgN):
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


def determine_patterns(pattern_set, patterns_to_equality, patterns_to_pgN):
    """
    Determine which patterns are useful in determining introgression
    Inputs:
    pattern_set -- a set containing all patterns of interest
    patterns_to_equality --- a mapping of site patterns to site patterns with equivalent p(gt|st)
    patterns_to_pgN --- a mapping of site patterns to their total p(g|N) value for a network
    Outputs:
    terms1 --- a set of patterns to count and add to each other to determine introgression
    terms2 --- a set of other patterns to count and add to each other to determine introgression
    """

    terms1 = set([])
    terms2 = set([])

    # Iterate over each pattern to determine the terms of interest
    for pattern1 in pattern_set:

        pat1_prob = patterns_to_pgN[pattern1]

        if pattern1 in patterns_to_equality.keys():
            for pattern2 in patterns_to_equality[pattern1]:

                pat2_prob = patterns_to_pgN[pattern2]

                if pat1_prob > pat2_prob:
                    terms1.add(pattern1)
                    terms2.add(pattern2)

                elif pat1_prob < pat2_prob:
                    terms1.add(pattern2)
                    terms2.add(pattern1)

    inverted1 = pattern_inverter(terms1)
    for pattern in inverted1:
        terms1.add(''.join(pattern))

    inverted2 = pattern_inverter(terms2)
    for pattern in inverted2:
        terms2.add(''.join(pattern))

    return terms1, terms2


def generate_statistic_string(patterns_of_interest):
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


##### Function for calculating statistic


def calculate_significance(left, right, verbose= False, alpha= 0.05):
    """
    Determines statistical significance based on a chi-squared goodness of fit test
    Input:
    left --- the total count for site patterns in the left term of the statistic
    right --- the total count for site patterns in the right term of the statistic
    Output:
    significant --- a boolean corresponding to whether or not the result is statistically significant
    """

    # Calculate the test statistic
    if left + right > 0:
        chisq = abs((left - right)**2 / float(left + right))
    else:
        chisq = 0

    # Calculate the p-value based on a chi square distribtion with df = 1
    pval = 1 - stats.chi2.cdf(chisq, 1)

    if pval < alpha:
        signif = True
    else:
        signif = False

    if verbose:
        return signif, left, right, chisq, pval
    else:
        return signif


def calculate_L(alignment, taxa_order, patterns_of_interest, verbose=False, alpha=0.05):
    """
    Calculates the L statistic for the given alignment
    Input:
    alignment --- a sequence alignment in phylip format
    taxa_order --- the desired order of the taxa
    patterns_of_interest --- a tuple containing the sets of patterns used for determining a statistic
    window_size --- the desired window size
    window_offset --- the desired offset between windows
    Output:
    l_stat --- the L statistic value
    significant --- a boolean denoting if the l_stat value is statistically significant
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

        # Statistic can only be calculated where the nucleotides are known
        if "-" not in bases and len(bases) == 2:

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
            site_string = pattern_string_generator([site_pattern])[0]

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

    # Verbose output
    if verbose:
        significant, left_counts, right_counts, chisq, pval = calculate_significance(terms1_total, terms2_total, verbose, alpha)
        return l_stat, significant, left_counts, right_counts, chisq, pval

    # Standard output
    else:
        significant = calculate_significance(terms1_total, terms2_total, alpha)
        return l_stat, significant


def calculate_windows_to_L(alignment, taxa_order, patterns_of_interest, window_size, window_offset, verbose= False, alpha=0.05):
    """
    Calculates the L statistic for the given alignment
    Input:
    alignment --- a sequence alignment in phylip format
    taxa_order --- the desired order of the taxa
    patterns_of_interest --- a tuple containing the sets of patterns used for determining a statistic
    window_size --- the desired window size
    windw_offset --- the desired offset between windows
    Output:
    l_stat --- the L statistic value
    windows_to_l --- a mapping of window indices to L statistic values
    """

    # Separate the patterns of interest into their two terms
    terms1 = patterns_of_interest[0]
    terms2 = patterns_of_interest[1]

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

    i = 0
    num_windows = 0
    if window_size > length_of_sequences:
        num_windows = 1
        window_size = length_of_sequences
    else:
        # Determine the total number of windows needed
        while (i + window_size - 1 < length_of_sequences):
            i += window_offset
            num_windows += 1

    site_idx = 0
    windows_to_l = {}

    # Iterate over each window
    for window in range(num_windows):

        terms1_counts = defaultdict(int)
        terms2_counts = defaultdict(int)

        # Iterate over the indices in each window
        for window_idx in range(window_size):

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

            # Statistic can only be calculated where the nucleotides are known
            if "-" not in bases and len(bases) == 2:

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
                site_string = pattern_string_generator([site_pattern])[0]

                # If the site string is a pattern of interest add to its count for one of the terms
                if site_string in terms1:
                    terms1_counts[site_string] += 1

                elif site_string in terms2:
                    terms2_counts[site_string] += 1

            # Increment the site index
            site_idx += 1

        terms1_total = sum(terms1_counts.values())
        terms2_total = sum(terms2_counts.values())

        numerator = terms1_total - terms2_total
        denominator = terms1_total + terms2_total

        if denominator != 0:
            l_stat = numerator / float(denominator)
        else:
            l_stat = 0

        # Verbose output
        if verbose:
            signif, left_counts, right_counts, chisq, pval = calculate_significance(terms1_total, terms2_total, verbose, alpha)
            # The line below can be changed to add more information to the windows to L mapping
            windows_to_l[window] = (l_stat, signif, chisq, pval)

        # Standard output
        else:
            signif = calculate_significance(terms1_total, terms2_total)
            windows_to_l[window] = (l_stat, signif)

        # Account for overlapping windows
        site_idx += (window_offset - window_size)

    return windows_to_l


##### Functions for total ordering


def branch_adjust(species_tree):
    """
    Create all possible combinations of branch lengths for the given species tree
    Input:
    species_tree --- a newick string containing the overall species tree
    Output:
    adjusted_trees --- a set of trees with all combinations of branch lengths
    """
    branch_lengths = [.5, 1.0, 2.0, 4.0]
    adjusted_trees = set([])

    taxa = []
    pattern = "((?<=\()[\w]+)|((?<=\,)[\w]+)"
    leaves = re.findall(pattern, species_tree)
    for leaf in leaves:
        if leaf[0] == '':
            taxa.append(leaf[1])
        else:
            taxa.append(leaf[0])

    for b in branch_lengths:
        new_t = species_tree
        for taxon in taxa:
            new_t = new_t.replace(taxon, "{0}:{1}".format(taxon, b))
        new_t = new_t.replace("),", "):{0},".format(b))
        adjusted_trees.add(new_t)

    return adjusted_trees, taxa


def network_branch_adjust(species_network):
    """
    Create all possible combinations of branch lengths for the given species network
    Input:
    species_tree --- a newick string containing the overall species network
    Output:
    adjusted_trees --- a set of trees with all combinations of branch lengths
    """

    branch_lengths = [.5, 1.0, 2.0, 4.0]
    adjusted_trees = set([])

    pattern = "((?<!\:)(\:\d+\.\d+))"
    lengths = re.findall(pattern, species_network)

    ############Adjust branch length stuff to account for all possible combinations

    for b in branch_lengths:
        new_t = species_network
        for l in lengths:
            new_t = new_t.replace(l[0], ":" + str(b))
        adjusted_trees.add(new_t)

    return list(adjusted_trees)


def network_adjust(species_network):
    """
    Create all possible combinations of inheritance probabilities for the given species network
    Input:
    species_network --- a newick string containing the overall species network
    Output:
    adjusted_networks --- a set of networks with all combinations of branch lengths
    """
    inheritance_probs = [0.1, 0.3]
    adjusted_networks = set([])

    pattern = "\:\:0\.\d+"
    reticulations = re.findall(pattern, species_network)

    for prob in inheritance_probs:
        new_net = species_network
        count = 0
        for r in reticulations:
            if count % 2 == 0:
                new_net = new_net.replace(r, "::{0}".format(str(prob)))
            else:
                new_net = new_net.replace(r, "::{0}".format(str(1 - prob)))
            count += 1
        net_set = network_branch_adjust(new_net)
        adjusted_networks = adjusted_networks.union(net_set)

    return adjusted_networks


def equality_sets(species_trees, network, taxa):
    """
    Create strings which represent the total ordering of p(gt|st)
    Input:
    species_tree --- a newick string containing the overall species tree without branch lengths
    Output:
    trees_to_equality --- a mapping of tree strings to a set of other trees with the same p(gt|st)
    trees_to_equality --- a mapping of tree strings to a set of other trees with the same p(gt|N)
    """
    st_to_pattern_probs = {}
    st_to_pattern_probs_N = {}
    trees_to_equality = {}
    trees_to_equality_N = {}

    outgroup = taxa[-1]
    gene_trees = generate_unique_trees(taxa, outgroup)

    newick_patterns = newicks_to_patterns_generator(taxa, gene_trees)

    for st in species_trees:
        ts_to_pgS, ts_to_pgN = calculate_newicks_to_stats(st, network, gene_trees)
        patterns_pgS, patterns_pgN = calculate_pattern_probabilities(newick_patterns, ts_to_pgS, ts_to_pgN)
        st_to_pattern_probs[st] = sorted(patterns_pgS.items(), key=lambda tup: tup[1], reverse=True)
        st_to_pattern_probs_N[st] = sorted(patterns_pgN.items(), key=lambda tup: tup[1], reverse=True)

    # Generate equality sets based on p(gt|st)
    for st in sorted(st_to_pattern_probs.keys()):

        gt_probs = st_to_pattern_probs[st]
        seen_trees = set([])

        for i in range(len(gt_probs)):

            gt1, prob1 = gt_probs[i]
            equal_trees = set([])

            seen = set([])
            for j in range(len(gt_probs)):

                gt2, prob2 = gt_probs[j]
                if prob1 == prob2 and gt1 != gt2 and gt2 not in seen_trees:
                    equal_trees.add(gt2)
                    seen.add(gt2)

            # Add the equality set to the mapping if tbe pattern is not already in the mapping and set is non empty
            if len(equal_trees) != 0 and gt1 not in seen_trees:
                if gt1 in trees_to_equality:
                    # If there is already an equality set for that pattern choose the larger of the two
                    if len(trees_to_equality[gt1]) > len(equal_trees):
                        break
                    else:
                        trees_to_equality[gt1] = equal_trees
                        seen_trees.add(gt1)
                        seen_trees = seen_trees.union(seen)
                else:
                    trees_to_equality[gt1] = equal_trees
                    seen_trees.add(gt1)
                    seen_trees = seen_trees.union(seen)

        # Generate equality sets based on p(gt|N)
        for st in sorted(st_to_pattern_probs_N.keys()):

            gt_probs = st_to_pattern_probs_N[st]
            seen_trees = set([])

            for i in range(len(gt_probs)):

                gt1, prob1 = gt_probs[i]
                equal_trees = set([])

                seen = set([])
                for j in range(len(gt_probs)):

                    gt2, prob2 = gt_probs[j]
                    if prob1 == prob2 and gt1 != gt2 and gt2 not in seen_trees:
                        equal_trees.add(gt2)
                        seen.add(gt2)

                # Add the equality set to the mapping if tbe pattern is not already in the mapping and set is non empty
                if len(equal_trees) != 0 and gt1 not in seen_trees:
                    if gt1 in trees_to_equality_N:
                        # If there is already an equality set for that pattern choose the larger of the two
                        if len(trees_to_equality_N[gt1]) > len(equal_trees):
                            break
                        else:
                            trees_to_equality_N[gt1] = equal_trees
                            seen_trees.add(gt1)
                            seen_trees = seen_trees.union(seen)
                    else:
                        trees_to_equality_N[gt1] = equal_trees
                        seen_trees.add(gt1)
                        seen_trees = seen_trees.union(seen)

    return trees_to_equality, trees_to_equality_N, patterns_pgS, patterns_pgN


def set_of_interest(trees_to_equality, trees_to_equality_N):
    """
    Inputs:
    trees_to_equality --- a mapping of tree strings to a set of other trees with the same p(gt|st)
    trees_to_equality_N --- a mapping of tree strings to a set of other trees with the same p(gt|N)
    Output:
    trees_of_interest --- a set of trees that changed equality under the species network
    """

    trees_of_interest = set([])

    for tree in trees_to_equality:

        if tree not in trees_to_equality_N:
            t_set = copy.deepcopy(trees_to_equality[tree])
            t_set.add(tree)
            trees_of_interest = trees_of_interest.union(t_set)
        elif trees_to_equality[tree] != trees_to_equality_N[tree]:
            t_set = copy.deepcopy(trees_to_equality[tree])
            t_set.add(tree)
            trees_of_interest = trees_of_interest.union(t_set)

    return trees_of_interest


def concat_directory(directory_path):
    """
        Concatenates all the alignments in a given directory and returns a single file.

        Input:
            directory_path --- a string path to the directory the use wants to use.
        Output:
            file_path --- a string path to the file that was created as a result of the concatenation.
    """

    # filter out hidden files
    filenames = filter(lambda n: not n.startswith(".") , natsorted(os.listdir(directory_path)))

    # get the number of lines on each file
    with open(os.path.join(directory_path, filenames[0]), "r") as f:
        n = len(list(f))

    # initialize a list with an empty string for each line
    output_file_list = [""] * n

    # Iterate over each folder in the given directory in numerical order
    for i in range(len(filenames)):

        # get full path of file
        input_file = os.path.join(directory_path, filenames[i])

        # if its a fasta file -> convert to phylip
        if filenames[i].endswith(".fa") or filenames[i].endswith(".fasta"):

            input_handle = open(input_file, "rU")
            output_handle = open(input_file + ".phylip", "w")

            alignments = AlignIO.parse(input_handle, "fasta")
            AlignIO.write(alignments, output_handle, "phylip-sequential")

            output_handle.close()
            input_handle.close()
            input_file = input_file + ".phylip"

        # create a list of the input files lines
        with open(input_file, 'r') as f:
            input_file_list = [l.rstrip() for l in list(f)]

        # print input_file_list

        for j in range(len(input_file_list)):
            # if this is the first file
            if i == 0:
                output_file_list[j] = input_file_list[j]
            else:
                if j == 0:
                    num_bp = int(input_file_list[0].split(" ")[2])
                    total_bp = int(output_file_list[j].split(" ")[2]) + num_bp
                    output_file_list[j] = " " + str(n - 1) + " " + str(total_bp)
                else:
                    output_file_list[j] += input_file_list[j].split(" ")[-1]


    # write the contents of the output file list to a text file
    with open(os.path.abspath(directory_path) + "/concatFile.phylip.txt", "w") as o:
        for line in output_file_list:
            print >> o, line

    return os.path.abspath(directory_path) + "/concatFile.phylip.txt"


def calculate_generalized(alignment, species_tree, reticulations, window_size, window_offset, verbose=False, alpha= 0.05, useDir=False, directory=""):
    """
    Calculates the L statistic for the given alignment
    Input:
        alignment --- a sequence alignment in phylip format
        taxa --- a list of the taxa in the desired order
        species_tree --- the inputted species tree over the given taxa
        reticulations --- a tuple containing two dictionaries mapping the start leaves to end leaves
        window_size --- the desired window size
        window_offset --- the desired offset between windows
        verbose --- a boolean for determining if extra information will be printed
        useDir --- a boolean for determining if the user wants to input an entire directory of alignments or only a single alignment
        directory --- a string path to the directory the use wants to use. NOTE: only necessary if useDir=True.
    Output:
        l_stat --- the L statistic value
    """

    print alignment

    st = re.sub("\:\d+\.\d+", "", species_tree)
    trees, taxa = branch_adjust(st)
    newick_patterns = newicks_to_patterns_generator(taxa, trees)
    network = generate_network_tree((0.1, 0.9), list(trees)[0], reticulations)
    trees_to_equality, trees_to_equality_N, patterns_pgS, patterns_pgN = equality_sets(trees, network, taxa)
    trees_of_interest = set_of_interest(trees_to_equality, trees_to_equality_N)
    increase, decrease = determine_patterns(trees_of_interest, trees_to_equality, patterns_pgN)

    if useDir:
        alignment = concat_directory(directory);

    if verbose:
        l_stat, significant, left_counts, right_counts, chisq, pval = calculate_L(alignment, taxa, (increase, decrease), verbose, alpha)
    else:
        l_stat, significant = calculate_L(alignment, taxa, (increase, decrease), verbose, alpha)

    windows_to_l = calculate_windows_to_L(alignment, taxa, (increase, decrease), window_size, window_offset, verbose, alpha)

    if verbose:
        print
        print "Newick strings with corresponding patterns: ", newick_patterns
        print
        print "Probability of gene tree patterns: ", patterns_pgS
        print
        print "Probability of species network patterns:", patterns_pgN
        print
        print "Patterns that were formerly equal with increasing probability: ", increase
        print "Patterns that were formerly equal with decreasing probability: ", decrease
        print
        print "Patterns of interest: ", increase, decrease
        print
        print "Statistic: ", generate_statistic_string((increase, decrease))
        print
        print "Overall Chi-Squared statistic: ", chisq
        print "Overall p value: ", pval
        print "Left term counts: ", left_counts
        print "Right term counts: ", right_counts
        print

    return l_stat, significant, windows_to_l


def plot_formatting(info_tuple, verbose=False):
    """
    Reformats and writes the dictionary output to a text file to make plotting it in Excel easy
    Input:
    info_tuple --- a triplet from the calculate_generalized output
    """

    l_stat, significance, windows_to_l = info_tuple

    num = 0
    file_name = "GeneralizedDResults{0}.txt".format(num)
    while os.path.exists(file_name):
        num += 1
        file_name ="GeneralizedDResults{0}.txt".format(num)

    with open(file_name, "w") as text_file:
        output_str = "Overall, {0}, {1} \n".format(l_stat, significance)
        text_file.write(output_str)
        for idx in windows_to_l:
            info = windows_to_l[idx]
            l_stat = info[0]
            significant = info[1]
            output_str = "{0}, {1}, {2} \n".format(idx, l_stat, significant)

            if verbose:
                chisq = info[2]
                pval = info[3]
                output_str = "{0}, {1}, {2}, {3}, {4} \n".format(idx, l_stat, significant, chisq, pval)

            text_file.write(output_str)



if __name__ == '__main__':
    # if we're running file directly and not importing it

    # Inputs for paper
    file = "C:\\Users\\travi\\Desktop\\concatFile.phylip.txt"
    species_tree = '((C,G),(((A,Q),L),R));'

    window_size, window_offset = 10000, 1000
    r = [('L', 'R')]
    plot_formatting(calculate_generalized(file, species_tree, r, window_size, window_offset, True))
    window_size, window_offset = 100000, 10000
    plot_formatting(calculate_generalized(file, species_tree, r, window_size, window_offset, True))

    window_size, window_offset = 10000, 1000
    r = [('Q', 'R')]
    plot_formatting(calculate_generalized(file, species_tree, r, window_size, window_offset, True))
    window_size, window_offset = 100000, 10000
    plot_formatting(calculate_generalized(file, species_tree, r, window_size, window_offset, True))

    window_size, window_offset = 10000, 1000
    r = [('Q', 'G')]
    plot_formatting(calculate_generalized(file, species_tree, r, window_size, window_offset, True))
    window_size, window_offset = 100000, 10000
    plot_formatting(calculate_generalized(file, species_tree, r, window_size, window_offset, True))

    # concat_directory("/Users/Peter/PycharmProjects/ALPHA/test_phylip_dir")
    # print calculate_generalized('/Users/Peter/PycharmProjects/ALPHA/CLFILE', '(((P1,P2),(P3,P4)),O);', [('P1', 'P3')], 50000, 50000, True)


        # file = 'C:\\Users\\travi\\Desktop\\clphylipseq.txt'
    # # r = [('L', 'R')]
    # r = [('Q', 'R')]
    # # r = [('Q', 'G')]
    # print calculate_generalized(file , '((C,G),(((A,Q),L),R));', r, 100000, 100000, True)

    # concat_directory("/Users/Peter/PycharmProjects/ALPHA/test_fasta_dir")
    # print calculate_generalized('/Users/Peter/PycharmProjects/ALPHA/CLFILE', '(((P1,P2),(P3,P4)),O);', [('P1', 'P3')], 50000, 50000, True)

    # species_tree, r = '(((P1:0.01,P2:0.01):0.01,(P3:0.01,P4:0.01):0.01):0.01,O:0.01);', [('P3', 'P1')]
    # species_tree = '(((P1,P2),(P3,P4)),O);'
    # alignment = "C:\\Users\\travi\\Documents\\PhyloVis\\exampleFiles\\ExampleDFOIL.phylip"
    # alignment = "C:\\Users\\travi\\Desktop\\seqfileNamed"
    # # lstat, signif, windows_to_l = calculate_generalized(alignment, species_tree, r, 1000, 1000, True, 0.05)
    # # plot_formatting((lstat, signif, windows_to_l))
    # plot_formatting(calculate_generalized('C:\\Users\\travi\\Desktop\\seqfileNamed', '(((P1,P2),(P3,P4)),O);', [('P3', 'P1')], 1000, 1000, True, 0.99), True)

    # print calculate_generalized('C:\\Users\\travi\\Desktop\\seqfileNamed', '(((P1,P2),(P3,P4)),O);', [('P1', 'P3')], 50000, 50000, True)

# python -c"from CalculateGeneralizedDStatistic import *; plot_formatting(calculate_generalized('C:\\Users\\travi\\Desktop\\seqfileNamed', '(((P1,P2),(P3,P4)),O);', [('P1', 'P3')], 100000, 100000, True, 0.01), True)"
