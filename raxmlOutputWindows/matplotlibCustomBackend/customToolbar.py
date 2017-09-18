import matplotlib
from customSubplotToolUi import CustomSubplotTool
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.qt_compat import QtCore, QtGui, QtWidgets
import os,six
import matplotlib.backends.qt_editor.figureoptions as figureoptions
import customFigureoptions

class CustomToolbar(NavigationToolbar):
    def __init__(self, canvas_, parent_):
        self.toolitems = (
            ('Home', 'Reset original view', 'home', 'home'),
            ('Back', 'Back to  previous view', 'back', 'back'),
            ('Forward', 'Forward to next view', 'forward', 'forward'),
            (None, None, None, None),
            ('Pan', 'Pan axes with left mouse, zoom with right', 'move', 'pan'),
            ('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom'),
            ('Subplots', 'Configure subplots', 'subplots', 'configure_subplots'),
            (None, None, None, None),
            ('Save', 'Save the figure', 'filesave', 'save_figure'),
        )

        NavigationToolbar.__init__(self, canvas_, parent_)

    def _icon(self, name):
        pm = QtGui.QPixmap(os.path.join(self.basedir, name))
        if hasattr(pm, 'setDevicePixelRatio'):
            pm.setDevicePixelRatio(self.canvas._dpi_ratio)
        return QtGui.QIcon(pm)

    def configure_subplots(self):
        image = os.path.join(matplotlib.rcParams['datapath'], 'images', 'matplotlib.png')
        dia = CustomSubplotTool(self.canvas.figure, self.parent)
        dia.setWindowIcon(QtGui.QIcon(image))
        dia.exec_()

    def generate_legend(self):
        allaxes = self.canvas.figure.get_axes()
        for ax in allaxes:
            leg = ax.legend()
            if leg:
                leg.draggable()
        self.canvas.draw()


    def edit_parameters(self):
        allaxes = self.canvas.figure.get_axes()
        if not allaxes:
            QtWidgets.QMessageBox.warning(self.parent, "Error", "There are no axes to edit.")
            return
        if len(allaxes) == 1:
            axes = allaxes[0]
        else:
            titles = []
            for axes in allaxes:
                name = (axes.get_title() or " - ".join(filter(None, [axes.get_xlabel(), axes.get_ylabel()])) or "<anonymous {} (id: {:#x})>".format(type(axes).__name__, id(axes)))
                titles.append(name)
            item, ok = QtWidgets.QInputDialog.getItem(self.parent, 'Customize', 'Select axes:', titles, 0, False)
            if ok:
                axes = allaxes[titles.index(six.text_type(item))]
            else:
                return

        customFigureoptions.figure_edit(axes, self)