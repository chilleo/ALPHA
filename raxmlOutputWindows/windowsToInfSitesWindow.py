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

class WindowsToInfSitesWindow(Window):
    def __init__(self, title, data, xLabel='Windows', yLabel='% Informative Sites'):
        Window.__init__(self, windowTitle='Windows to Informative Sites', legend=False)

        self.plotter.stat_scatter(data, title, xLabel=xLabel, yLabel=yLabel)
        self.show()


if __name__ == '__main__': # only runs if not imported

    # create a new instance of QApplication
    app = QtGui.QApplication(sys.argv)

    a = {0:0, 2:2, 4:3}

    # create window and plot
    form = WindowsToInfSitesWindow(a)

    # execute the app
    sys.exit(app.exec_())
