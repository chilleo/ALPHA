from standardWindow import Window
from PyQt4 import QtGui
import sys

"""
MS Percent Matching Window
~
Chabrielle Allen
Travis Benedict
Peter Dulworth
"""

class MSPercentMatchingWindow(Window):
    def __init__(self, title1, data1, xLabel1='', yLabel1='', groupLabels1=()):
        Window.__init__(self, windowTitle='Percent Sites Matching')

        self.plotter.barPlot(title1, data1, xLabel=xLabel1, yLabel=yLabel1, groupLabels=groupLabels1, subplotPosition=111)

        self.show()


if __name__ == '__main__': # only runs if not imported

    # create a new instance of QApplication
    app = QtGui.QApplication(sys.argv)

    a = {0:0, 2:2, 4:3}

    # create window and plot
    form = MSPercentMatchingWindow(a)

    # execute the app
    sys.exit(app.exec_())
