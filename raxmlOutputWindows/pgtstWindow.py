from standardWindow import Window
from PyQt4 import QtGui
import sys


"""
P(gene tree | species tree)
~
Chabrielle Allen
Travis Benedict
Peter Dulworth
"""


class PGTSTWindow(Window):
    def __init__(self, windowsToPGTST, title, xLabel='', yLabel=''):
        Window.__init__(self, windowTitle='P(gene tree | species tree)', legend=False)

        self.plotter.stat_scatter(windowsToPGTST, title, xLabel, yLabel)
        self.show()


if __name__ == '__main__': # only runs if not imported

    # create a new instance of QApplication
    app = QtGui.QApplication(sys.argv)

    a = {0:0, 1:0, 2:0.5, 3:0.6, 4:0.3, 5:0.5}

    # create window and plot
    form = PGTSTWindow(a, "p(gt|st)", "Windows", "Probability")
    form.show()

    # execute the app
    sys.exit(app.exec_())
