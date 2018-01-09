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

class LStatisticWindow(Window):
    def __init__(self, lWindows, title='Windows to L Statistic', xLabel='Window Indices', yLabel='L Statistic values'):
        Window.__init__(self, windowTitle='L Statistic Window')

        self.plotter.stat_scatter(lWindows, title, xLabel=xLabel, yLabel=yLabel)
        self.show()


if __name__ == '__main__': # only runs if not imported

    # create a new instance of QApplication
    app = QtGui.QApplication(sys.argv)

    a = {0:0, 2:2, 4:3}

    # create window and plot
    form = LStatisticWindow(a)

    # execute the app
    sys.exit(app.exec_())
