import sys
from module import plotter
import matplotlib
matplotlib.use('Qt4Agg')  # necessary for mac pls don't remove -- needs to be before pyplot is imported but after matplotlib is imported
from matplotlib import pyplot as plt
from matplotlibCustomBackend.customToolbar import CustomToolbar
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.qt_compat import QtCore, QtGui


class Window(QtGui.QMainWindow):
    def __init__(self, windowTitle='Window', x=0, y=0, legend=True, parent=None):
        super(Window, self).__init__(parent)

        # set UI style -- options: u'Windows', u'Motif', u'CDE', u'Plastique', u'Cleanlooks', u'Macintosh (aqua)'
        # QtGui.QApplication.setStyle(QtGui.QStyleFactory.create(u'Plastique'))

        # layout
        self.setAutoFillBackground(True)
        self.centralwidget = QtGui.QWidget(self)
        self.verticalLayout = QtGui.QVBoxLayout(self.centralwidget)
        self.verticalLayout.setMargin(0)
        self.setCentralWidget(self.centralwidget)
        self.setBackgroundColor(QtCore.Qt.white)

        # create menu bar
        self.menubar = QtGui.QMenuBar(self)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 340, 22))

        # create menus
        self.menuFile = QtGui.QMenu(self.menubar)
        self.menuFile.setTitle("File")
        self.menuConfigurePlot = QtGui.QMenu(self.menubar)
        self.menuConfigurePlot.setTitle("Configure Plot")

        # create actions
        self.actionSaveAs = QtGui.QAction(self)
        self.actionSaveAs.setText("Save As...")
        self.actionConfigureSubplots = QtGui.QAction(self)
        self.actionConfigureSubplots.setText('Configure Subplots')
        self.actionConfigureAxis = QtGui.QAction(self)
        self.actionConfigureAxis.setText('Configure Axis and Curve')

        # add actions to menu
        self.menuFile.addAction(self.actionSaveAs)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menuConfigurePlot.addAction(self.actionConfigureSubplots)
        self.menuConfigurePlot.addAction(self.actionConfigureAxis)
        self.menubar.addAction(self.menuConfigurePlot.menuAction())

        # enable menu bar
        self.setMenuBar(self.menubar)

        QtCore.QMetaObject.connectSlotsByName(self)

        # get arguments
        self.windowTitle = windowTitle

        # moves menu bar into application -- mac only windows sux
        self.menubar.setNativeMenuBar(False)

        # set window title
        self.setWindowTitle(self.windowTitle)

        self.initCanvas()

        self.toolbar = CustomToolbar(self.canvas, self)
        self.verticalLayout.addWidget(self.toolbar)

        self.setWindowPosition(x, y)

        # bind file menu actions
        self.connect(self.actionSaveAs, QtCore.SIGNAL('triggered()'), self.toolbar.save_figure)

        # bind configure menu actions
        self.connect(self.actionConfigureSubplots, QtCore.SIGNAL('triggered()'), self.toolbar.configure_subplots)
        self.connect(self.actionConfigureAxis, QtCore.SIGNAL('triggered()'), self.toolbar.edit_parameters)

        if legend:
            self.generateLegendMenu()

        # create instance of Plotter class
        self.plotter = plotter.Plotter()

    def initCanvas(self):
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.verticalLayout.addWidget(self.canvas)

    def generateLegendMenu(self):
        self.menuLegend = QtGui.QMenu(self.menubar)
        self.menuLegend.setTitle('Legend')
        self.actionGenerateLegend = QtGui.QAction(self)
        self.actionGenerateLegend.setText('Generate Draggable Legend')
        self.menuLegend.addAction(self.actionGenerateLegend)
        self.menubar.addAction(self.menuLegend.menuAction())
        self.connect(self.actionGenerateLegend, QtCore.SIGNAL('triggered()'), self.generateLegend)

    def generateLegend(self):
        self.toolbar.generate_legend()

    def setBackgroundColor(self, color):
        """
            change background color to white
        """

        p = self.palette()
        p.setColor(self.backgroundRole(), color)
        self.setPalette(p)

    def setWindowPosition(self, x, y):
        """
            positions the window relative to the top left corner of the screen (px)
        """
        self.move(x, y)

    def setWindowSize(self, x, y):
        """
            sets size of window
        """
        self.resize(x, y)

    def closeEvent(self, QCloseEvent):
        # plt.clf()
        self.emit(QtCore.SIGNAL("WINDOW_CLOSED"))
        print self.windowTitle + ' Closed'

    # def moveEvent(self, QMoveEvent):
    #     print self.fileName, self.pos()


if __name__ == '__main__':
    """
        code is executed if file is run directly -- i.e. not imported
    """

    # A new instance of QApplication
    app = QtGui.QApplication(sys.argv)

    # initialize main input window
    form = Window(windowTitle='Standard Window')
    form.show()

    # form.plotter.barPlot('bar', [1,2,3,4], groupLabels=[1,2,3,4])

    form.setWindowSize(600, 600)

    # and execute the app
    sys.exit(app.exec_())
