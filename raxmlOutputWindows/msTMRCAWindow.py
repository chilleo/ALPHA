from standardWindow import Window
from PyQt4 import QtGui
import sys

"""
MS TMRCA Window
~
Chabrielle Allen
Travis Benedict
Peter Dulworth
"""

class MSTMRCAWindow(Window):
    def __init__(self, data, labels):
        Window.__init__(self, windowTitle='TMRCA Graph')

        self.plotter.tmrca_graph(data, labels)
        self.show()


if __name__ == '__main__': # only runs if not imported

    # create a new instance of QApplication
    app = QtGui.QApplication(sys.argv)

    a = {0:0, 2:2, 4:3}
    b = [1,2,3]

    # create window and plot
    form = MSTMRCAWindow(a, b)

    # execute the app
    sys.exit(app.exec_())
