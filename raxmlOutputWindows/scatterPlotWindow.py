import sys
from standardWindow import Window
from PyQt4 import QtGui

"""
    Scatter Plot Window
    ~
    Chabrielle Allen
    Travis Benedict
    Peter Dulworth
"""


class ScatterPlotWindow(Window):
    def __init__(self, title, windowsToTopologies, colors, heights):
        Window.__init__(self, windowTitle='Windows to Top Topologies')

        # plot
        self.plotter.topologyScatter(title, windowsToTopologies, colors, heights)
        self.show()


if __name__ == '__main__': # only runs if not imported

    # create a new instance of QApplication
    app = QtGui.QApplication(sys.argv)

    a = {0: '(C,(G,O),H);', 1: '((C,G),O,H);', 2: '(C,(G,O),H);', 3: '(C,(G,O),H);', 4: '(C,(G,O),H);', 5: '(C,(G,O),H);', 6: '(C,(G,O),H);', 7: '(C,(G,O),H);', 8: '((C,G),O,H);', 9: '(C,(G,O),H);'}
    b = ['#ff0000', '#0000ff', '#ff0000', '#ff0000', '#ff0000', '#ff0000', '#ff0000', '#ff0000', '#0000ff', '#ff0000']
    c = [1, 0, 1, 1, 1, 1, 1, 1, 0, 1]

    # create window and plot
    form = ScatterPlotWindow('', a, b, c)
    form.show()
    form.plot()

    # execute the app
    sys.exit(app.exec_())

