from standardWindow import Window
from PyQt4 import QtGui
import sys

"""
Robinson Foulds Window
~
Chabrielle Allen
Travis Benedict
Peter Dulworth
"""

class RobinsonFouldsWindow(Window):
    def __init__(self, title1, data1, title2=None, data2=None, xLabel='Windows', yLabel='RF Distance'):
        Window.__init__(self, windowTitle='Robinson Foulds Distance From MS Truth', legend=False)

        if title2 == None:
            self.plotter.stat_scatter(data1, title1, xLabel, yLabel, subplotPosition=111)
        else:
            self.plotter.stat_scatter(data1, title1, xLabel, yLabel, subplotPosition=211)
            self.plotter.stat_scatter(data2, title2, xLabel, yLabel, subplotPosition=212)

        self.show()


if __name__ == '__main__': # only runs if not imported

    # create a new instance of QApplication
    app = QtGui.QApplication(sys.argv)

    a = {0:0, 2:2, 4:3}

    # create window and plot
    form = RobinsonFouldsWindow(a)

    # execute the app
    sys.exit(app.exec_())
