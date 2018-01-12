from matplotlib import pyplot as plt
import numpy as np
import math

def calculate_d(alignment, window_size, window_offset, taxon1, taxon2, taxon3, taxon4):
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

    windows_to_d = {}

    sequence_list = []
    taxon_list = []

    with open(alignment) as f:

        # Create a list of each line in the file
        lines = f.readlines()

        # First line contains the number and length of the sequences
        first_line = lines[0].split()
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
            i += window_offset
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
        if denominator_window != 0:
            d_window = numerator_window / float(denominator_window)
        else:
            d_window = 0

        # Map the window index to its D statistic
        windows_to_d[window] = d_window

        # Reset numerator and denominator
        numerator_window = 0
        denominator_window = 0

        # Account for overlapping windows
        site_idx += (window_offset - window_size)

    d_numerator = 0
    d_denominator = 0

    # Iterate over the site indices
    for site_idx in range(length_of_sequences):

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

        d_numerator += (ABBA - BABA)
        d_denominator += (ABBA + BABA)

    if d_denominator != 0:
        d_stat = d_numerator / float(d_denominator)
    else:
        d_stat = 0

    return d_stat, windows_to_d


def stat_scatter(dataMap, dstat, filename, title='Windows to D Statistic', xLabel='Window Indices', yLabel='D Statistic values', subplotPosition=111):
    """
        Creates a scatter plot with the x-axis being the
        windows and the y-axis being the statistic to
        be graphed.

        Input:
            i. dataMap -- a mapping
            ii. dstat -- the overall d statistic calculation
            iii. filename -- the filename for the outputted plot
            iv. name -- the name of the save file
            v. title -- the title of the plot
            vi. xLabel -- the label for the x axis
            vii. yLabel -- the label for the y axis
            
    """

    ax = plt.subplot(subplotPosition)
    ax.set_title(title, fontsize=15)

    # sizes plot circles
    circleArea = math.pi * (3) ** 2

    # makes x values integers
    xValues = dataMap.keys()
    xIntegerValues = []

    for x in xValues:
        xIntegerValues.append(int(x))

    x = np.array(xIntegerValues)

    # get y values from dictionary and convert to an np array
    yValues = dataMap.values()
    y = np.array(yValues)

    # create scatter plot
    ax.scatter(x, y, s=circleArea, c='#000000', alpha=1)
    ax.axhline(y=dstat, color='r', linestyle='-')

    ax.set_ylim((-1, 1))

    # label the axes
    ax.set_xlabel(xLabel, fontsize=10)
    ax.set_ylabel(yLabel, fontsize=10)

    plt.savefig(filename)


def plot_d(alignment, window_size, window_offset, taxon1, taxon2, taxon3, taxon4, filename):
    """
    
    """

    dstat, windows_to_d = calculate_d(alignment, window_size, window_offset, taxon1, taxon2, taxon3, taxon4)
    min_d = min(windows_to_d.values())
    max_d = max(windows_to_d.values())
    range_d = abs(max_d - min_d)
    stat_scatter(windows_to_d, dstat, filename)
    return range_d


# species_tree, r = '(((P1:0.01,P2:0.01):0.01,(P3:0.01,P4:0.01):0.01):0.01,O:0.01);', [('P3', 'P1'),('P3', 'P2')]
alignment = "C:\\Users\\travi\\Documents\\PhyloVis\\exampleFiles\\Example.phylip"
filename = "C:\\Users\\travi\\Desktop\\test.png"
# print plot_d(alignment, 50000, 50000, "P1", "P2", "P3", "O", filename)

# python -c "from CalculateDStatistic import *; print calculate_d('AlignmentPhylipFile', WindowSize, WIndowOffset,'taxon1', 'taxon2', 'taxon3', 'outgroup')"
# python -c "from CalculateDStatistic import *; print plot_d('C:\\Users\\travi\\Documents\\PhyloVis\\exampleFiles\\Example.phylip', 50000, 50000,'P1', 'P2', 'P3', 'O')"

