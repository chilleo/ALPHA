from collections import defaultdict
from natsort import natsorted
import os
import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from PyQt4 import QtCore

"""
Functions:
    __init__(self, parent=None)
    is_site_informative(self, site)    
    calculate_informativeness(self, window_directory, window_offset)
    line_graph_generator(self, dictionary, xlabel, ylabel, name)
    heat_map_generator(self, dictionary, name)
~
Chabrielle Allen
Travis Benedict
Peter Dulworth
"""


class InformativeSites(QtCore.QThread):
    def __init__(self, parent=None):
        super(InformativeSites, self).__init__(parent)

    def is_site_informative(self, site):
        """
        Determines if a site is informative or not
        Input:
        site --- a list of bases located at a site in the alignment
        Output:
        1 if a site is informative 0 if a site is uninformative
        """

        # Create a mapping of bases to the number of times they occur
        base_to_counts = defaultdict(int)

        # Iterate over each base in the list site
        for base in site:

            # Add one each time a base occurs
            base_to_counts[base] += 1

        # Create a list of counts in descending order
        base_counts = sorted(base_to_counts.values(), reverse=True)

        if len(base_counts) >= 2:
            # If two different bases occur at least twice the site is informative
            if (base_counts[0] >= 2) and (base_counts[1] >= 2):
                return 1

            else:
                return 0
        else:
            return 0

    def calculate_informativeness(self, window_directory, window_offset, percentage):
        """
        Calculates information about informative sites in an alignment
        Input:
        window_directory --- the location of the folder containing the phylip window files
        window_offset --- the offset that was used to create the windows
        percentage --- the percent of the total alignment to look at the site index is scaled accordingly
        Output:
        sites_to_informative --- a mapping of each site in the alignment to 1 if informative 0 if not
        windows_to_informative_count --- a mapping of each window number to the number of informative sites it has
        windows_to_informative_pct --- a mapping of each window number to the percentage of informative sites it has
        pct_informative --- the percentage of informative sites over the entire alignment
        """

        # Initialize the site index to 0
        site_idx = 0

        # Represent percentage as a decimal
        pct = float(percentage) / 100

        sites_to_informative = defaultdict(int)
        windows_to_informative_count = defaultdict(int)
        windows_to_informative_pct = {}
        total_window_size = 0

        # Iterate over each folder in the given directory in numerical order
        for filename in natsorted(os.listdir(window_directory)):

            # If file is a phylip file get the number of the window
            if filename.endswith(".phylip"):
                file_number = filename.replace("window", "")
                file_number = int(file_number.replace(".phylip", ""))

                input_file = os.path.join(window_directory, filename)

                sequence_list = []

                with open(input_file) as f:

                    # Create a list of each line in the file
                    lines = f.readlines()

                    # First line contains the number and length of the sequences
                    first_line = lines[0].split()
                    number_of_sequences = int(first_line[0])
                    length_of_sequences = int(first_line[1])

                for line in lines[1:]:
                    # Add each sequence to a list
                    sequence = line.split()[1]
                    sequence_list.append(sequence)

                # Increment based on the percentage of the alignment desired
                increment = int(math.ceil((length_of_sequences * pct) / length_of_sequences))

                # Iterate over the indices in each window
                for window_idx in range(length_of_sequences):

                    site = []

                    # Iterate over each sequence in the alignment
                    for sequence in sequence_list:

                        # Add each base in a site to a list
                        site.append(sequence[window_idx])

                    # Determine if a site is informative
                    informative = self.is_site_informative(site)

                    # If the site has not been visited before add to mappings (deals with overlapping windows)
                    if site_idx not in sites_to_informative and informative:
                        # If the site is informative add 1 to the mappings otherwise add 0
                        sites_to_informative[site_idx] += informative

                    windows_to_informative_count[file_number] += informative

                    # Increment the site index
                    site_idx += increment

                # Account for overlapping windows
                site_idx += (window_offset - length_of_sequences)

                # Map windows_to_informative_count to a percentage
                windows_to_informative_pct[file_number] = windows_to_informative_count[file_number] *\
                                                          (100 / float(length_of_sequences))

                total_window_size += length_of_sequences

        total_num_informative = sum(windows_to_informative_count.values())

        # Add in the last site index if it is not already in the informative mapping
        # This is useful for plotting reasons
        if site_idx not in sites_to_informative:
            sites_to_informative[site_idx] = 0

        if total_window_size == 0:
            pct_informative = 0
        else:
            pct_informative = float(total_num_informative * 100) / total_window_size

        return sites_to_informative, windows_to_informative_count, windows_to_informative_pct, pct_informative


if __name__ == '__main__':  # if we're running file directly and not importing it
    # travys window dir
    # window_dir = "C:\\Users\\travi\\Documents\\Evolutionary-Diversity-Visualization-Python\\windows"

    # peters window dir
    # window_dir = '/Users/Peter/PycharmProjects/Evolutionary-Diversity-Visualization-Python/windows'

    # chabs window dir ?
    # window_dir = ''

    infs = InformativeSites()

    sites_to_informative, windows_to_informative_count, windows_to_informative_pct, pct_informative = infs.calculate_informativeness(window_dir, 50000)

    # print str(pct_informative) + "%"

    infs.line_graph_generator(windows_to_informative_pct, "Windows", "Percentage of Informative Sites", "pctInformative.png")

    infs.heat_map_generator(sites_to_informative, "HeatMapInfSites.png")
