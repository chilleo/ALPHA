from standardWindow import Window
from PyQt4 import QtGui
import sys

"""
All Trees Window
~
Chabrielle Allen
Travis Benedict
Peter Dulworth
"""


class AllTreesWindow(Window):
    def __init__(self, title, colorScheme, rooted=False, outGroup=False):
        Window.__init__(self, windowTitle='All Trees Window', legend=False)

        self.plotter.topologyColorizer(title, colorScheme, rooted=rooted, outgroup=outGroup)
        self.show()


if __name__ == '__main__': # only runs if not imported

    # create a new instance of QApplication
    app = QtGui.QApplication(sys.argv)

    a = {'((C,G),O,H);': '#0000ff', '(C,(G,O),H);': '#ff0000', '((C,G),(O,H));': '#00ff00'}

    # create window and plot
    form = AllTreesWindow('', a, rooted=False, outGroup=False)
    form.show()
    form.plot()

    # execute the app
    sys.exit(app.exec_())
