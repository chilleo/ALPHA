import matplotlib
matplotlib.use('Qt4Agg')  # necessary for mac pls don't remove -- needs to be before pyplot is imported but after matplotlib is imported
from matplotlib import pyplot as plt
from PyQt4 import QtCore
import math, re, numpy as np
from Bio import Phylo
from ete3 import Tree, TreeNode
from cStringIO import StringIO
from Bio.Graphics.GenomeDiagram import GraphSet
from Bio.Graphics import GenomeDiagram
from reportlab.lib import colors


class Plotter(QtCore.QThread):
    def __init__(self, parent=None):
        super(Plotter, self).__init__(parent)

        # list of colors for plots
        self.COLORS = ['#ff0000', '#0000ff', '#32cd32', '#ba55d3', '#87cefa', '#ffa500', '#ff1493', '#a020f0', '#00ced1', '#adff2f', '#ffd700', '#1e90ff', '#ff7f50', '#008000', '#ffc0cb', '#8a2be2']

        # list of patterns for line styles
        self.PATTERNS = ['-', '--', ':', '-.']

        self.plot = 'genomeAtlas'

    def topologyDonut(self, title, labels, sizes, donut_colors, subplotPosition=111):
        """
            Creates a donut chart showing the breakdown of the top 'num' topologies.

            Inputs:
            i. labels -- a list of labels outputted by top_freqs()[1]
            ii. sizes  -- a list of sizes outputted by top_freqs()[2]
            iii. donut_colors -- a list of colors outputted by donut_colors()

            Returns:
                A donut chart with the number of times a topology occurs and 'Other Topologies' for topologies that occur less than the
                most frequent 'num' topologies as the labels, and a list tops of the top 'num' scores.
        """

        ax = plt.subplot(subplotPosition, aspect='equal')
        ax.set_title(title, fontsize=15)

        ax.pie(sizes, explode=None, labels=labels, colors=donut_colors, autopct=None, shadow=False)

        # create circle behind pie chart to outline it
        outer_circle = plt.Circle((0, 0), 1, color='#000000', fill=False, linewidth=1.25)

        # impose circle over pie chart to make a donut chart
        inner_circle = plt.Circle((0, 0), 0.65, color='#000000', fc='#ffffff', linewidth=1.25)

        ax.add_artist(inner_circle)
        ax.add_artist(outer_circle)

        return ax

    def topologyScatter(self, title, wins_to_tops, scatter_colors, y, subplotPosition=111):
        """
            Creates a scatter plot showing the topology as the y-axis and the window as the x-axis.

            Input:
                i. wins_to_tops   -- window to topology mapping outputted by windows_to_newick()[0]
                ii. scatter_colors -- list of colors outputted by topology_colors()[1]
                iii. y          -- list of y-axis values outputted by topology_colors()[2]

            Returns:
                A scatter plot with topologies as the x-axis and windows as the y-axis.
        """

        ax = plt.subplot(subplotPosition)
        ax.set_title(title, fontsize=15)

        # area of plotted circles
        circleArea = math.pi * (3) ** 2

        # size y-axis on plot
        ax.set_yticks(np.arange(len(wins_to_tops) + 1, 0))

        # x-axis is window numbers
        windows = wins_to_tops.keys()

        # for each index, and unique topology
        for (i, topology) in enumerate(set(wins_to_tops.values())):
            xc = [top for (j, top) in enumerate(windows) if wins_to_tops.values()[j] == topology]
            yc = [top for (j, top) in enumerate(y) if wins_to_tops.values()[j] == topology]
            cols = [c for (j, c) in enumerate(scatter_colors) if wins_to_tops.values()[j] == topology]
            ax.scatter(xc, yc, s=circleArea, c=cols, alpha=1, edgecolors='#000000', label=topology)
            ax.grid = True

        # labels axes
        ax.set_xlabel('Windows', fontsize=10)
        ax.set_ylabel('Top Newick Strings', fontsize=10)

        return ax

    def stat_scatter(self, dataMap, title, xLabel, yLabel, subplotPosition=111):
        """
            Creates a scatter plot with the x-axis being the
            windows and the y-axis being the statistic to
            be graphed.

            Input:
                i. dataMap -- a mapping
                ii. name -- the name of the save file
                iii. title -- the title of the plot
                iv. xLabel -- the label for the x axis
                v. yLabel -- the label for the y axis

            Returns:
                A scatter plot with windows as the x-axis and a statistic as the y-axis.
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

        # label the axes
        ax.set_xlabel(xLabel, fontsize=10)
        ax.set_ylabel(yLabel, fontsize=10)

        return ax

    def barPlot(self, title, data, xLabel='', yLabel='', groupLabels=(), heightLabels=False, subplotPosition=111):
        """
            generates bar chart
        """

        ax = plt.subplot(subplotPosition)

        ax.set_title(title, fontsize=15)
        ax.set_xlabel(xLabel, fontsize=10)
        ax.set_ylabel(yLabel, fontsize=10)
        ax.set_xticks([])

        numberOfBars = len(data)
        indices = range(numberOfBars)  # the x locations for the groups
        width = .667  # the width of the bars
        bars = []

        if heightLabels:
            ax.set_yticks([])
            for i in indices:
                bars.append(ax.bar(i, data[i], width, label=groupLabels[i]))
                for bar in bars[i]:
                    h = bar.get_height()
                    ax.text(bar.get_x() + bar.get_width() / 2.0, h, '%d' % int(h), ha='center', va='bottom')
        else:
            for i in indices:
                ax.bar(i, data[i], width, label=groupLabels[i])

        ax.set_xmargin(0.1)
        # ax.legend()

        if len(groupLabels) > 0:
            ax.set_xticks(indices)
            ax.set_xticklabels(groupLabels)

        return ax

    def heatMap(self, title, mapping, subplotPosition=111):
        """
            Create a heat map based on the inputted dictionary

            Input:
                i. title --- a string for the title of the plot
                ii. mapping --- a dictionary mapping integers to floats or integers
        """

        ax = plt.subplot(subplotPosition)
        ax.set_title(title, fontsize=15)

        y = mapping.values()
        N = len(y)
        x = range(N)
        width = 30 * (N / float(20000))

        plt.bar(x, y, width, color="black")

        ax.spines['top'].set_position(("data", 1)) #Reposition top spine
        ax.spines['left'].set_bounds(0, 1) #Shorten left and right spines
        ax.spines['right'].set_bounds(0, 1)
        cur_axes = plt.gca()
        cur_axes.axes.get_yaxis().set_visible(False)

        self.emit(QtCore.SIGNAL('HEATMAP_COMPLETE'))

        return ax

    def treeImage(self, title, newick, rooted=False, outgroup=False):
        """
            Given a newick string, creates an image of the tree.
            Used in L Statistic GUI.
        """

        plt.title(title, fontsize=15)
        plt.axis('off')

        ax = plt.subplot(1, 1, 1)
        ax.text(0, 0, str(1) + ' Occurrences', style='italic')

        # disable axis
        ax.axis('off')

        # Create the tree object
        tree = Phylo.read(StringIO(newick), "newick")
        tree.rooted = rooted

        if rooted:
            tree.root_with_outgroup(outgroup)

        # Create the tree image
        Phylo.draw(tree, axes=ax, do_show=False)

    def topologyColorizer(self, title, newicksToColors, topologies_to_counts, rooted=False, outgroup=False):
        """
            Create colored tree topology images based on a color scheme where the color of a tree is determined by the frequency that it occurs.

            Inputs:
                i. color scheme --- a dictionary mapping newick strings to colors
                ii. rooted --- a boolean parameter corresponding to the tree being rooted
                iii. outgroup --- a string of the desired taxon to root at
        """

        plt.title(title, fontsize=15)
        plt.axis('off')

        # count number of top topologies
        numTopTopologies = 0
        for newick in newicksToColors:
            if newick != "Other":
                numTopTopologies += 1

        # create a count for the number of the topologies
        count = 1
        # Iterate over each newick string in color_scheme
        for newick in newicksToColors:
            newick_count = topologies_to_counts[newick]
            if newick != "Other":

                if numTopTopologies < 4:
                    ax = plt.subplot(numTopTopologies, 1, count)
                elif numTopTopologies == 5:
                    order = [None, 1,3,5,7,9]
                    ax = plt.subplot(3, 3, order[count])
                else:
                    ax = plt.subplot(int(round(numTopTopologies / 2.0)), 2, count)

                ax.text(0, 0, str(newick_count) + ' Occurences' , style='italic')

                # disable axis
                ax.axis('off')

                # Create the tree object and assign it to the appropriate color
                tree = Phylo.read(StringIO(newick), "newick")
                tree.rooted = rooted

                if rooted:
                    tree.root_with_outgroup(outgroup)

                tree.root.color = newicksToColors[newick]

                # Create the tree image
                Phylo.draw(tree, axes=ax, do_show=False)

                count += 1

    def doubleLineGraph(self, list1, list2, confidenceThreshold, xLabel='', yLabel=''):
        """
            Create a line graph based on the inputted dictionary

            Input:
                i. list1 --- a list of integers
                ii. list2 --- a list of integers of equal length to list1
                iii. xLabel --- a string for the labeling the x-axis
                iv. yLabel --- a string for the labeling the y-axis
                v. name --- a string for the image name

        """

        ax = plt.subplot(111)
        ax.set_title('Number of Internal Nodes After Contraction Confidence Threshold: ' + str(confidenceThreshold))

        rangeX = range(len(list1))

        ax.plot(rangeX, list1, "-", label='Before Contraction')
        ax.plot(rangeX, list2, "-", label='After Contraction')

        ax.set_xlabel(xLabel)
        ax.set_ylabel(yLabel)

    def nlineGraph(self, dataSet, title, xLabel='', yLabel='', subplotPosition=111, groupLabels1=()):
        """
            Create a line graph based on the inputted dictionary

            Input:
                i. dataSet --- a list of dictionaries mapping integers to floats or integers
                ii. xLabel --- a string for the labeling the x-axis
                iii. yLabel --- a string for the labeling the y-axis
                iv. title --- title of plot
                v. subplotPosition - matplotlib position specifying dimensions of subplot grid and position of subplot
        """

        ax = plt.subplot(subplotPosition)

        ax.set_title(title, fontsize=15)
        ax.set_xlabel(xLabel)
        ax.set_ylabel(yLabel)

        for data in dataSet:
            x = data.keys()
            y = data.values()

            ax.plot(x, y, "-")

        ax.legend(groupLabels1)

        return ax

    def lineGraph(self, data, title, xLabel='', yLabel='', subplotPosition=111):
        """
            Create a line graph based on the inputted dictionary

            Input:
                i. data --- a dictionary mapping integers to floats or integers
                ii. xLabel --- a string for the labeling the x-axis
                iii. yLabel --- a string for the labeling the y-axis
                iv. title --- title of plot
                v. subplotPosition - matplotlib position specifying dimensions of subplot grid and position of subplot
        """

        ax = plt.subplot(subplotPosition)

        ax.set_title(title, fontsize=15)
        ax.set_xlabel(xLabel)
        ax.set_ylabel(yLabel)

        x = data.keys()
        y = data.values()

        ax.plot(x, y, "-")

        return ax

    def tmrca_graph(self, sites_to_newick_mappings, labels, topology_only=False, subplotPosition=111):
        """
            Plots a line graph comparing tree heights from different MS files.

            Inputs:
                i. sites_to_newick_mappings -- a list of the mappings outputted by sites_to_newick_ms()
                ii. topology_only: If set to True, distance between nodes will be referred to the number of nodes between them.
                    In other words, topological distance will be used instead of branch length distances.

            Returns:
                i. A line graph with the tree height as the y-axis and the site number as the x-axis.
        """

        print labels

        ax = plt.subplot(subplotPosition)

        ax.set_title('TMRCA Line Graph')
        ax.set_xlabel('SNP Site Number')
        ax.set_ylabel('TMRCA')

        trees = []
        roots = []
        leaves = []
        dist = []
        heights = []

        # iterate over each mapping in list
        for i in range(len(sites_to_newick_mappings)):
            mapping = sites_to_newick_mappings[i]
            for tree in mapping:
                # iterate over mapping to get trees
                trees.append(mapping[tree])

            for j in range(len(trees)):
                # get tree roots
                roots.append(Tree.get_tree_root(Tree(trees[j])))

                # get distance from roots to farthest leaves
                leaves.append(TreeNode.get_farthest_leaf(roots[j], topology_only))

            for k in range(len(leaves)):
                # regular expression to get height values from list of farthest leaves
                dist.append(re.findall(', \d{1,}.\d{1,}', str(leaves[k])))

                # format with regular expression to remove unnecessary tokens
                heights.append(re.sub("\[', |']", '', str(dist[k])))

            # resets ind to prevent index error in linestyle pattern
            # if i > 3:
            #     ind = random.randint(0, 3)
            # else:
            #     ind = i

            # plot line graph
            ax.plot(sites_to_newick_mappings[0].keys(), heights, c=self.COLORS[i], linestyle=self.PATTERNS[i % len(self.PATTERNS)], label=labels[i])

            # clear lists
            trees = []
            roots = []
            leaves = []
            dist = []
            heights = []

        leg = ax.legend()
        if leg:
            leg.draggable()

        return ax

    def generateCircleGraph(self, file, windows_to_top_topologies, topologies_to_colors, window_size, window_offset, sites_to_informative):
        """
            Creates genetic circle graph showing the windows and the areas where each topology appears

            Inputs:
                i. file --- phylip file inputted in GUI
                ii. windows_to_top_topologies --- mapping outputted by windows_to_newick()[0]
                iii. topologies_to_colors --- mapping outputted by topology_colors()[0]
                iv. window_size --- size inputted in GUI
                v. window_offset --- size inputted in GUI
            Returns:
                A genetic circle graph GenomeAtlas.
        """

        ############################# Format Data #############################

        # get the length of the sequence
        with open(file) as f:
            length_of_sequences = int(f.readline().split()[1])
        f.close()

        # accounts for offset
        windows_to_top_topologies2 = {}
        for window in windows_to_top_topologies:
            windows_to_top_topologies2[window * window_offset] = windows_to_top_topologies[window]
        windows_to_top_topologies = windows_to_top_topologies2.items()

        # gets windows
        windows = []
        for window_topology in windows_to_top_topologies:
            windows.append(window_topology[0])

        for i in range(length_of_sequences):
            if i not in windows:
                windows_to_top_topologies.append((i, 0))

        # maps data to topologies
        topologies_to_data = {}
        for topology in topologies_to_colors:
            topologies_to_data[topology] = []

        for topology in topologies_to_colors:
            for window in windows_to_top_topologies:
                if topology == window[1]:
                    topologies_to_data[topology].append(tuple([window[0], 1]))
                else:
                    topologies_to_data[topology].append(tuple([window[0], -0.1]))

        # maps data to colors
        data_to_colors = {}
        for topology in topologies_to_data:
            data_to_colors[str(topologies_to_data[topology])] = topologies_to_colors[topology]

        includeOther = False
        if 'Other' in topologies_to_data:
            includeOther = True
            # -1 because top_topologies_to_colors includes 'Other'
            number_of_top_topologies = len(topologies_to_colors) - 1
        else:
            number_of_top_topologies = len(topologies_to_colors)

        # removes 'Other' from mapping
        if includeOther:
            minor_topology_data = topologies_to_data['Other']
            del topologies_to_data['Other']
        data = topologies_to_data.values()

        # separates data into windowed and un-windowed
        windowed_data = []
        nonwindowed_data = []
        for i in range(length_of_sequences):
            if i not in windows:
                nonwindowed_data.append(tuple([i, 1]))
                windowed_data.append(tuple([i, 0]))
            else:
                for j in range(i, i + window_size):
                    windowed_data.append((j, 1))
                nonwindowed_data.append(tuple([i, 0]))

        ############################# Build Graph #############################

        # name of the figure
        name = "../plots/GenomeAtlas"
        graphStyle = 'bar'

        # create the diagram -- highest level container for everything
        diagram = GenomeDiagram.Diagram(name)

        diagram.new_track(
            1,
            greytrack=0,
            name="Track",
            height=2,
            hide=0,
            scale=1,
            scale_color=colors.black,
            scale_font='Helvetica',
            scale_fontsize=6,
            scale_fontangle=45,
            scale_ticks=1,
            scale_largeticks=0.3,
            scale_smallticks=0.1,
            scale_largetick_interval=(length_of_sequences / 6),
            scale_smalltick_interval=(length_of_sequences / 12),
            scale_largetick_labels=1,
            scale_smalltick_labels=0
        )

        if includeOther:
            diagram \
                .new_track(2, name="Minor Topologies", height=1.0, hide=0, greytrack=0, greytrack_labels=2,
                           greytrack_font_size=8, grey_track_font_color=colors.black, scale=1, scale_ticks=0,
                           axis_labels=0) \
                .new_set('graph') \
                .new_graph(minor_topology_data, style=graphStyle,
                           colour=colors.HexColor(data_to_colors[str(minor_topology_data)]),
                           altcolour=colors.transparent, linewidth=1)

        for i in range(number_of_top_topologies):
            # create tracks -- and add them to the diagram
            if i == 0:
                diagram \
                    .new_track(i + 3, name="Track" + str(i + 1), height=1.0, hide=0, greytrack=0,
                               greytrack_labels=2, greytrack_font_size=8, grey_track_font_color=colors.black,
                               scale=1, scale_ticks=0, axis_labels=0) \
                    .new_set('graph') \
                    .new_graph(data[i], style=graphStyle, colour=colors.HexColor(data_to_colors[str(data[i])]),
                               altcolour=colors.transparent, linewidth=1)
            else:
                diagram \
                    .new_track(i + 3, name="Track" + str(i + 1), height=1.0, hide=0, greytrack=0,
                               greytrack_labels=2, greytrack_font_size=8, grey_track_font_color=colors.black,
                               scale=0) \
                    .new_set('graph') \
                    .new_graph(data[i], style=graphStyle, colour=colors.HexColor(data_to_colors[str(data[i])]),
                               altcolour=colors.transparent, linewidth=1)

        # outer ring shit
        graph_set = GraphSet('graph')
        graph_set.new_graph(nonwindowed_data, style=graphStyle, color=colors.HexColor('#cccccc'),
                            altcolour=colors.transparent)
        graph_set.new_graph(windowed_data, style=graphStyle, color=colors.HexColor('#2f377c'),
                            altcolour=colors.transparent)

        # diagram \
        #     .new_track(i + 5, name="Track" + str(i + 2), height=2, hide=0, greytrack=0, greytrack_labels=2,
        #                greytrack_font_size=8, grey_track_font_color=colors.black, scale=0) \
        #     .add_set(graph_set)

        ############################################

        diagram \
            .new_track(i + 4, name="Track" + str(i + 1), height=1.0, hide=0, greytrack=0,
                       greytrack_labels=2, greytrack_font_size=8, grey_track_font_color=colors.black,
                       scale=0) \
            .new_set('graph') \
            .new_graph(sites_to_informative.items(), style='line', colour=colors.coral,
                       altcolour=colors.transparent, linewidth=0.01)

        ###########################################################

        diagram.draw(format="circular", pagesize='A5', orientation='landscape', x=0.0, y=0.0, track_size=1.88,
                     tracklines=0, circular=0, circle_core=0.3, start=0, end=length_of_sequences - 1)

        # save the file
        # diagram.write("../plots/GenomeAtlas.eps", "EPS")
        # diagram.write("plots/GenomeAtlas.svg", "SVG")
        diagram.write("plots/GenomeAtlas.pdf", "PDF")
        diagram.write("plots/GenomeAtlas.png", "PNG")

        self.emit(QtCore.SIGNAL('CIRCLE_GRAPH_COMPLETE'))

    def run(self):
        if self.plot == 'genomeAtlas':
            self.generateCircleGraph(self.alignment, self.windowsToTopTopologies, self.topologiesToColors, self.windowSize, self.windowOffset, self.sitesToInformative)
        elif self.plot == 'heatmap':
            self.heatMap(self.title, self.sitesToInformative)


if __name__ == '__main__':  # if we're running file directly and not importing it
    p = Plotter()
    # p.figureBarPlot([1,2,3,4], 'name')
    # a = {0: '(C,(G,O),H);', 1: '((C,G),O,H);', 2: '(C,(G,O),H);', 3: '(C,(G,O),H);', 4: '(C,(G,O),H);', 5: '(C,(G,O),H);', 6: '(C,(G,O),H);', 7: '(C,(G,O),H);', 8: '((C,G),O,H);', 9: '(C,(G,O),H);'}
    # b = ['#ff0000', '#0000ff', '#ff0000', '#ff0000', '#ff0000', '#ff0000', '#ff0000', '#ff0000', '#0000ff', '#ff0000']
    # c = [1, 0, 1, 1, 1, 1, 1, 1, 0, 1]
    # p.heatMap('title', )


    # p.barPlot('', [1,2,3,4,5 ], '', '', groupLabels=('one', 'two', '3', 4, '5'))
    p.treeImage("title", "(a,(b,c))")

    plt.show()