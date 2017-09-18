from standardWindow import Window
from PyQt4 import QtGui
import sys

"""
MS Robinson Foulds Window
~
Chabrielle Allen
Travis Benedict
Peter Dulworth
"""

class MSRobinsonFouldsWindow(Window):
    def __init__(self, title, data, xLabel1='', yLabel1='', groupLabels1=()):
        Window.__init__(self, windowTitle='Robinson Foulds Distance From MS Truth')

        self.plotter.nlineGraph(data, title, xLabel=xLabel1, yLabel=yLabel1, subplotPosition=111)

        self.show()


if __name__ == '__main__': # only runs if not imported

    # create a new instance of QApplication
    app = QtGui.QApplication(sys.argv)

    a = {0:0, 2:2, 4:3}

    # create window and plot
    form = MSRobinsonFouldsWindow(a)

    # execute the app
    sys.exit(app.exec_())
