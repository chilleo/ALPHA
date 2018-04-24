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
    def __init__(self, lWindows, sigArray):
        Window.__init__(self, windowTitle='L Statistic Window')

        self.plotter.sorted_scatter(lWindows, sigArray, 'Windows to Generalized D', 'Windows', 'Generalized D Values')
        self.show()


if __name__ == '__main__': # only runs if not imported

    # create a new instance of QApplication
    app = QtGui.QApplication(sys.argv)

    a = {0:0, 2:2, 4:3}
    b = [0, 0, 1]

    # create window and plot
    form = LStatisticWindow(a, b)

    # execute the app
    sys.exit(app.exec_())
