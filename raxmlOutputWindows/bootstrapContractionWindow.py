from standardWindow import Window
from PyQt4 import QtGui
import sys

"""
Bootstrap Contraction Window
~
Chabrielle Allen
Travis Benedict
Peter Dulworth
"""


class BootstrapContractionWindow(Window):
    def __init__(self, seq1, seq2, confidenceThreshold, xLabel='', yLabel=''):
        Window.__init__(self, windowTitle='Bootstrap Contraction Window')

        self.plotter.doubleLineGraph(seq1, seq2, confidenceThreshold, xLabel=xLabel, yLabel=yLabel)
        self.show()


if __name__ == '__main__': # only runs if not imported

    # create a new instance of QApplication
    app = QtGui.QApplication(sys.argv)

    # create window and plot
    form = BootstrapContractionWindow([0,1,2,3,4],[3,2,3,2,3], 70, xLabel="Window Indices", yLabel="Number of Internal Nodes")
    form.show()

    # execute the app
    sys.exit(app.exec_())
