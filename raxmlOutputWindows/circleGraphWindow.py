from standardWindow import Window
from PyQt4 import QtGui, QtCore
import sys
import matplotlib
matplotlib.use('Qt4Agg')  # necessary for mac pls don't remove -- needs to be before pyplot is imported but after matplotlib is imported
from matplotlib import pyplot as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.image as mpimg
import svgutils.compose as sc


"""
Functions:
    __init__(self)
~
Chabrielle Allen
Travis Benedict
Peter Dulworth
"""


class CircleGraphWindow(Window):
    def __init__(self, alignment, windowsToTopTopologies, topologiesToColors, windowSize, windowOffset, sitesToInformative):
        Window.__init__(self, windowTitle='Genome Atlas Window')

        self.connect(self.plotter, QtCore.SIGNAL('CIRCLE_GRAPH_COMPLETE'), self.paintCircleGraph)

        self.plotter.plot = 'genomeAtlas'
        self.plotter.alignment = alignment
        self.plotter.windowsToTopTopologies = windowsToTopTopologies
        self.plotter.topologiesToColors = topologiesToColors
        self.plotter.windowSize = windowSize
        self.plotter.windowOffset = windowOffset
        self.plotter.sitesToInformative = sitesToInformative

        self.plotter.start()

    def initCanvas(self):
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.verticalLayout.addWidget(self.canvas)
        self.ax = plt.subplot(111)
        self.ax.axis('off')

    def paintCircleGraph(self):
        img = mpimg.imread('plots/GenomeAtlas.png')
        self.ax.imshow(img)
        self.show()

    def generateLegend(self):
        """
            todo: implement a custom function to generate a legend for circle graph
        """
        pass

if __name__ == '__main__': # only runs if not imported

    # create a new instance of QApplication
    app = QtGui.QApplication(sys.argv)

    alignment = '../testFiles/phylip.txt'
    windowSize = 10
    windowOffset = 10

    # create window and plot
    form = CircleGraphWindow()

    # execute the app
    sys.exit(app.exec_())
