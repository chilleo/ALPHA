from standardWindow import Window
from PyQt4 import QtGui
import sys

"""
D-Statistic Window
~
Chabrielle Allen
Travis Benedict
Peter Dulworth
"""

class DStatisticWindow(Window):
    def __init__(self, dWindows, title='Windows to D-Statistic', xLabel='Window Indices', yLabel='D-Statistic values'):
        Window.__init__(self, windowTitle='D-Statistic Window')

        self.plotter.stat_scatter(dWindows, title, xLabel=xLabel, yLabel=yLabel)
        self.show()


if __name__ == '__main__': # only runs if not imported

    # create a new instance of QApplication
    app = QtGui.QApplication(sys.argv)

    a = {0:0, 2:2, 4:3}

    # create window and plot
    form = DStatisticWindow(a)

    # execute the app
    sys.exit(app.exec_())
