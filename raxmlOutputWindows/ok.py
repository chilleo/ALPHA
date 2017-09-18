import sys
from PyQt4 import QtGui, QtCore
import matplotlib
matplotlib.use("qt4agg")
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

class GUI(QtGui.QMainWindow):
    def __init__(self, parent=None):
        super(GUI, self).__init__(parent)
        self.buildLayout()
        self.buildMenus()
        self.menuBar()
        self.statusBar()

        ## Style Sheets
        self.splitter.setStyleSheet("QSplitter::handle:horizontal {background-color:   #ccc}")
        self.controlWidget.setStyleSheet(".QWidget {background-color:   #0ff}")
        menuStyle = """.QMenuBar {background-color:   #0ff}
            QMenuBar::item {background: transparent}
            QMenuBar::item:selected {background: #8ff}"""
        self.statusBar().setStyleSheet(".QStatusBar {background-color:   #0ff}")
        self.menuBar().setStyleSheet(menuStyle)
        # .....THIS DOESN"T WORK !! .....
        self.mplFig.setStyleSheet("QWidget {background-color:   #f00}")
        self.mplFig.setStyleSheet("QWidget {background:   #f00}")
        self.mplFig.setStyleSheet("QWidget {color:   #f00}")
        # self.mplFig.setStyleSheet("""QWidget {
        #          background-color: #0f0;
        #         }
        #
        #      QWidget::item {
        #          background: #0f0;
        #      }""")


    def buildLayout(self):
        self.controlWidget = QtGui.QWidget(self)
        self.plotList  = QtGui.QListWidget(self)
        self.combo  = QtGui.QComboBox(self)
        self.button = QtGui.QPushButton('Plot')
        self.combo.addItems(['1','2','3','4'])
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.plotList)
        layout.addWidget(self.combo)
        layout.addWidget(self.button)
        self.controlWidget.setLayout(layout)
        self.mplFig  = MplGrapher()
        self.splitter = QtGui.QSplitter(QtCore.Qt.Horizontal)
        self.splitter.addWidget(self.controlWidget)
        self.splitter.addWidget(self.mplFig)
        self.setCentralWidget(self.splitter)
        # QtGui.QApplication.setStyle(QtGui.QStyleFactory.create('Plastique'))
    def buildMenus(self):
        openFile = QtGui.QAction('Open', self)
        self.fileMenu = self.menuBar().addMenu('&File')
        self.fileMenu.addAction(openFile)

class MplGrapher(QtGui.QWidget):
    def __init__(self,parent=None):
        super(MplGrapher, self).__init__(parent)
        self.initFigure()

    def initFigure(self):
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.navbar = NavigationToolbar(self.canvas, self)
        self.figure.add_subplot(1,1,1)
        self.layout = QtGui.QVBoxLayout()
        self.layout.addWidget(self.navbar)
        self.layout.addWidget(self.canvas)
        self.setLayout(self.layout)

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    main = GUI()
    main.show()
    sys.exit(app.exec_())