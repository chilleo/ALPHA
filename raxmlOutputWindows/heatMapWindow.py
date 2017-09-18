from standardWindow import Window
from PyQt4 import QtGui, QtCore
import sys
import matplotlib
matplotlib.use('Qt4Agg')  # necessary for mac pls don't remove -- needs to be before pyplot is imported but after matplotlib is imported
from matplotlib import pyplot as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas

"""
Informative Sites Heatmap
~
Chabrielle Allen
Travis Benedict
Peter Dulworth
"""


class HeatMapWindow(Window):
    def __init__(self, title, sitesToInformative):
        Window.__init__(self, windowTitle='Informative Sites Heatmap', legend=False)

        self.connect(self.plotter, QtCore.SIGNAL('HEATMAP_COMPLETE'), self.show)

        self.plotter.plot = 'heatmap'
        self.plotter.title = title
        self.plotter.sitesToInformative = sitesToInformative

        self.plotter.start()

    def initCanvas(self):
        self.figure = plt.figure(figsize=(15, 2))
        self.canvas = FigureCanvas(self.figure)
        self.verticalLayout.addWidget(self.canvas)

if __name__ == '__main__': # only runs if not imported

    # create a new instance of QApplication
    app = QtGui.QApplication(sys.argv)

    a = {0:0, 2:2, 4:3}

    # create window and plot
    form = HeatMapWindow('Heat Map', a)

    # execute the app
    sys.exit(app.exec_())
