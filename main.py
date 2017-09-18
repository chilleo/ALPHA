# utilities
import sip, sys, os, re, webbrowser
sip.setapi('QString', 2)
from PyQt4 import QtGui, QtCore
from functools import partial

# GUI
from raxmlOutputWindows import allTreesWindow, donutPlotWindow, scatterPlotWindow, circleGraphWindow, pgtstWindow, robinsonFouldsWindow, heatMapWindow, bootstrapContractionWindow, dStatisticWindow, msRobinsonFouldsWindow, msPercentMatchingWindow, msTMRCAWindow, windowsToInfSitesWindow
from module import gui_layout as gui

# logic
from module import RAxMLOperations as ro
from module import topologyPlots as tp
from module import statisticCalculations as sc
from module import fileConverterController as fc
from module import informativeSites as infSites
from module import bootstrapContraction as bc
from module import msComparison as ms


class PhyloVisApp(QtGui.QMainWindow, gui.Ui_PhylogeneticVisualization):
    def __init__(self, parent=None):
        super(PhyloVisApp, self).__init__(parent)

        # remove any leftover files from previous raxml trials
        badFileNames = ['RAxML_result', 'RAxML_randomTree', 'RAxML_log', 'RAxML_info', 'RAxML_bestTree', 'RAxML_bipartitions', 'RAxML_bipartitionsBranchLabels', 'RAxML_bootstrap']
        for fileName in os.listdir('.'):
            nameWithoutExtension = os.path.splitext(fileName)[0]
            for file in badFileNames:
                if nameWithoutExtension == file:
                    os.remove(fileName)

        # if 'plots' folder doesn't exist -> create it
        if not os.path.isdir('plots'):
            os.mkdir('plots')

        # remove all files in plots folder
        for fileName in os.listdir('plots'):
            os.remove('plots/' + fileName)

        # initialize gui_layout
        self.setupUi(self)

        # set UI style -- options: u'Windows', u'Motif', u'CDE', u'Plastique', u'Cleanlooks', u'Macintosh (aqua)'
        QtGui.QApplication.setStyle(QtGui.QStyleFactory.create(u'Macintosh (aqua)'))

        self.dStatisticTaxonComboBoxes = [self.dTaxonComboBox1, self.dTaxonComboBox2, self.dTaxonComboBox3, self.dTaxonComboBox4]
        self.raxmlTaxonComboBoxes = [self.outgroupComboBox]
        self.speciesTreeComboBoxes = [self.speciesTreeComboBox]

        # moves menu bar into application -- mac only windows sux
        self.menubar.setNativeMenuBar(False)

        # set GUI icon
        self.setWindowIcon(QtGui.QIcon('imgs/phylovisLogo.png'))

        # self.welcomeLogoImage.setScaledContents(True)
        self.welcomeLogoImage.setPixmap(QtGui.QPixmap('imgs/phylovisLogo.png'))

        # create new instance of RaxmlOperations class
        self.raxmlOperations = ro.RAxMLOperations()
        # create new instance of TopologyPlotter class
        self.topologyPlotter = tp.TopologyPlotter()
        # create new instance of Statistics Calculations class
        self.statisticsCalculations = sc.StatisticsCalculations()
        # create new instance of Informative Sites class
        self.informativeSites = infSites.InformativeSites()
        # create new instance of BootstrapContraction class
        self.bootstrapContraction = bc.BootstrapContraction()
        # create new instance of MsComparison class
        self.msComparison = ms.MsComparison()
        # create new instance of FileConverter class
        self.fileConverter = fc.FileConverter()

        # mapping from: windows --> page index
        self.windows = {'welcomePage': 0, 'inputPageRax': 1, 'inputPageFileConverter': 2, 'inputPageMS': 3, 'inputPageDStatistic': 4}
        # mapping from: windows --> dictionary of page dimensions
        self.windowSizes = {'welcomePage': {'x': 459, 'y': 245}, 'inputPageRax': {'x': 600, 'y': 600}, 'inputPageFileConverter': {'x': 459, 'y': 403}, 'inputPageMS': {'x': 600, 'y': 680}, 'inputPageDStatistic': {'x': 600, 'y': 600}}
        # mapping from: windows --> dictionary of page dimensions
        self.windowLocations = {'welcomePage': {'x': 600, 'y': 300}, 'inputPageRax': {'x': 500, 'y': 175}, 'inputPageFileConverter': {'x': 600, 'y': 300}, 'inputPageMS': {'x': 520, 'y': 100}, 'inputPageDStatistic': {'x': 500, 'y': 175}}
        # mapping from: mode --> page
        self.comboboxModes_to_windowNames = {'RAx_ML': 'inputPageRax', 'File Converter': 'inputPageFileConverter', 'MS Comparison': 'inputPageMS', 'D Statistic': 'inputPageDStatistic'}
        # mapping from: mode --> menu action
        self.comboboxModes_to_actionModes = {'RAx_ML': self.actionRax, 'File Converter': self.actionConverter, 'MS Comparison': self.actionMS, 'D Statistic': self.actionDStatistic}
        # if running windows os
        if sys.platform == 'win32':
            # mapping from: windows --> dictionary of page dimensions
            self.windowSizes = {'welcomePage': {'x': 459, 'y': 245}, 'inputPageRax': {'x': 925, 'y': 688}, 'inputPageFileConverter': {'x': 630, 'y': 375}, 'inputPageMS': {'x': 675, 'y': 815}, 'inputPageDStatistic': {'x': 600, 'y': 570}}

        # default values
        self.runComplete = False
        self.checkboxWeighted.setEnabled(False)
        self.outgroupComboBox.setEnabled(False)
        self.outgroupLabel.setEnabled(False)
        self.bootstrapGroupBox.setEnabled(False)
        self.outgroupGroupBox.setEnabled(False)
        self.speciesTreeOutGroupGroupBox.setEnabled(False)
        self.dStatisticLabel.setEnabled(False)
        self.speciesTreeRaxmlCommandEntry.setEnabled(False)
        self.customRaxmlCommandEntry.setEnabled(False)
        self.progressBar.reset()
        self.generateSpeciesTreeProgressBar.reset()
        self.rooted = False
        self.stackedWidget.setCurrentIndex(0)
        self.raxmlToolBox.setCurrentIndex(0)
        self.raxmlOptionsTabWidget.setCurrentIndex(0)
        self.resize(self.windowSizes['welcomePage']['x'], self.windowSizes['welcomePage']['y'])
        self.outputFileConverterEntry.setText(os.getcwd())

        # open documentation
        self.actionDocumentation.triggered.connect(lambda: self.openURL('https://peterdulworth.github.io/PhyloVis'))

        # only allow integers in the following fields
        self.setValidator(self.windowSizeEntry, 'Int')
        self.setValidator(self.windowOffsetEntry, 'Int')
        self.setValidator(self.numberOfTopTopologiesEntry, 'Int')
        self.setValidator(self.confidenceLevelEntry, 'Int')
        self.setValidator(self.numberOfBootstrapsEntry, 'Int')
        self.setValidator(self.msWindowSizeEntry, 'Int')
        self.setValidator(self.msWindowOffsetEntry, 'Int')
        self.setValidator(self.dWindowSizeEntry, 'Int')
        self.setValidator(self.dWindowOffsetEntry, 'Int')

        # **************************** RAXML PAGE ****************************#

        # selecting a mode in the menu bar -> deselects all other modes first
        # change the input mode based on which mode is selected in the menu bar
        self.actionRax.triggered.connect(lambda: self.ensureSingleModeSelected(self.actionRax, 'inputPageRax'))
        self.actionConverter.triggered.connect(lambda: self.ensureSingleModeSelected(self.actionConverter, 'inputPageFileConverter'))
        self.actionMS.triggered.connect(lambda: self.ensureSingleModeSelected(self.actionMS, 'inputPageMS'))
        self.actionDStatistic.triggered.connect(lambda: self.ensureSingleModeSelected(self.actionDStatistic, 'inputPageDStatistic'))

        # triggers select file dialogs
        self.inputFileBtn.clicked.connect(lambda: self.getFileName(self.inputFileEntry))
        self.newickFileBtn.clicked.connect(lambda: self.getFileName(self.newickFileEntry))

        # regenerates each graph every time checkbox is checked
        self.checkboxCircleGraph.stateChanged.connect(lambda: self.updatedDisplayWindows(btnClicked=self.checkboxCircleGraph))
        self.checkboxScatterPlot.stateChanged.connect(lambda: self.updatedDisplayWindows(btnClicked=self.checkboxScatterPlot))
        self.checkboxAllTrees.stateChanged.connect(lambda: self.updatedDisplayWindows(btnClicked=self.checkboxAllTrees))
        self.checkboxDonutPlot.stateChanged.connect(lambda: self.updatedDisplayWindows(btnClicked=self.checkboxDonutPlot))
        self.checkboxHeatMap.stateChanged.connect(lambda: self.updatedDisplayWindows(btnClicked=self.checkboxHeatMap))
        self.checkboxWindowsToInfSites.stateChanged.connect(lambda: self.updatedDisplayWindows(btnClicked=self.checkboxWindowsToInfSites))
        self.checkboxPGTST.stateChanged.connect(lambda: self.updatedDisplayWindows(btnClicked=self.checkboxPGTST))
        self.checkboxRobinsonFoulds.stateChanged.connect(lambda: self.updatedDisplayWindows(btnClicked=self.checkboxRobinsonFoulds))

        # toggle what inputs are actionable based on checkboxes
        self.checkboxRobinsonFoulds.clicked.connect(lambda: self.toggleEnabled(self.checkboxWeighted))
        self.checkboxRooted.stateChanged.connect(lambda: self.toggleEnabled(self.outgroupComboBox))
        self.checkboxRooted.stateChanged.connect(lambda: self.toggleEnabled(self.outgroupLabel))
        self.checkboxBootstrap.stateChanged.connect(lambda: self.toggleEnabled(self.bootstrapGroupBox))
        self.checkboxRooted.stateChanged.connect(lambda: self.toggleEnabled(self.outgroupGroupBox))
        self.checkBoxCustomRaxml.stateChanged.connect(lambda: self.toggleEnabled(self.customRaxmlCommandEntry))
        self.checkboxSpeciesTreeRooted.stateChanged.connect(lambda: self.toggleEnabled(self.speciesTreeOutGroupGroupBox))
        self.checkboxSpeciesTreeUseCustomRax.stateChanged.connect(lambda: self.toggleEnabled(self.speciesTreeRaxmlCommandEntry))

        # RAxML Events
        self.connect(self.inputFileEntry, QtCore.SIGNAL('FILE_SELECTED'), lambda: self.updateTaxonComboBoxes(self.raxmlTaxonComboBoxes, self.inputFileEntry))
        self.connect(self.inputFileEntry, QtCore.SIGNAL('FILE_SELECTED'), lambda: self.updateTaxonComboBoxes(self.speciesTreeComboBoxes, self.inputFileEntry))
        self.connect(self.raxmlOperations, QtCore.SIGNAL('RAX_PER'), self.progressBar.setValue)
        self.connect(self.raxmlOperations, QtCore.SIGNAL('RAX_COMPLETE'), self.raxmlComplete)
        self.connect(self.raxmlOperations, QtCore.SIGNAL('RAX_COMPLETE'), self.updatedDisplayWindows)
        self.connect(self.raxmlOperations, QtCore.SIGNAL('SPECIES_TREE_PER'), self.generateSpeciesTreeProgressBar.setValue)
        self.connect(self.raxmlOperations, QtCore.SIGNAL('SPECIES_TREE_COMPLETE'), partial(self.message, type='Err'))
        self.connect(self.raxmlOperations, QtCore.SIGNAL('INVALID_ALIGNMENT_FILE'), lambda: self.message('Invalid File', 'Invalid alignment file. Please choose another.', 'Make sure your file has 4 sequences and is in the phylip-relaxed format.', type='Err'))

        # self.connect(self.topologyPlotter, QtCore.SIGNAL('CIRCLE_GRAPH_COMPLETE'), lambda: self.openWindow(self.circleGraphWindow))

        # run RAX_ML and generate graphs
        self.runBtn.clicked.connect(self.runRAxML)
        self.generateSpeciesTreeBtn.clicked.connect(self.generateSpeciesTree)

        # **************************** WELCOME PAGE ****************************#

        self.launchBtn.clicked.connect(self.initializeMode)

        # **************************** CONVERTER PAGE ****************************#

        self.fileTypeDocumentationBtn.clicked.connect(lambda: self.openURL('http://biopython.org/wiki/AlignIO'))

        self.fileConverterBtn.clicked.connect(lambda: self.getFileName(self.fileConverterEntry))
        self.outputFileConverterBtn.clicked.connect(lambda: self.openDirectory(self.outputFileConverterEntry))
        self.runFileConverterBtn.clicked.connect(lambda: self.convertFile())

        self.connect(self.fileConverter, QtCore.SIGNAL('FILE_CONVERTER_COMPLETE'), lambda: self.fileConverterProgressBar.setValue(100))
        self.connect(self.fileConverter, QtCore.SIGNAL('FILE_CONVERTER_COMPLETE'), self.message)
        self.connect(self.fileConverter, QtCore.SIGNAL('FILE_CONVERTER_ERR'), self.message)

        # **************************** MS PAGE ****************************#

        self.msCompareBtn.clicked.connect(self.runMSCompare)
        self.msFileBtn.clicked.connect(lambda: self.getFileName(self.msFileEntry))
        self.msSecondFileBtn.clicked.connect(lambda: self.getFileName(self.msSecondFileEntry))

        self.connect(self.msComparison, QtCore.SIGNAL('MS_COMPLETE'), self.plotMSCompare)
        self.connect(self.msComparison, QtCore.SIGNAL('MS_PER'), self.msProgressBar.setValue)
        self.connect(self.msComparison, QtCore.SIGNAL('MS_ERR'), self.message)

        self.checkboxCompareAgainstMS.clicked.connect(lambda: self.toggleEnabled(self.msMSCompareGroupBox))
        self.checkboxCompareAgainstRaxml.clicked.connect(lambda: self.toggleEnabled(self.msRaxmlCompareGroupBox))

        self.msRaxmlDirectoryBtn.clicked.connect(lambda: self.openDirectory(self.msRaxmlDirectoryEntry))

        self.msUploadAnother.clicked.connect(lambda: self.addFileEntry('msAdditionalFileHorizontalLayout', 'msAdditionalFileEntry', 'msAdditionalFileBtn', 'msRemoveFileBtn'))

        # **************************** D STATISTIC PAGE ****************************#

        # set background image
        self.imagePixmap = QtGui.QPixmap('imgs/tree.png')
        self.imageLabel.setScaledContents(True)
        self.imageLabel.setPixmap(self.imagePixmap)

        # run
        self.dAlignmentBtn.clicked.connect(lambda: self.getFileName(self.dAlignmentEntry))

        # when file entry text is changed
        self.connect(self.dAlignmentEntry, QtCore.SIGNAL("FILE_SELECTED"), lambda: self.updateTaxonComboBoxes(self.dStatisticTaxonComboBoxes, self.dAlignmentEntry, require4Taxons=True))

        # update progress bar
        self.connect(self.statisticsCalculations, QtCore.SIGNAL('D_PER'), self.dProgressBar.setValue)
        self.connect(self.statisticsCalculations, QtCore.SIGNAL('D_FINISHED'), self.displayDStatistic)

        # run
        self.dRunBtn.clicked.connect(self.runDStatistic)
        self.connect(self.statisticsCalculations, QtCore.SIGNAL('INVALID_ALIGNMENT_FILE'), partial(self.message, type='Err'))

    # **************************** WELCOME PAGE ****************************#

    def initializeMode(self):
        self.ensureSingleModeSelected(self.comboboxModes_to_actionModes[self.modeComboBox.currentText()], self.comboboxModes_to_windowNames[self.modeComboBox.currentText()])

    # **************************** D STATISTIC PAGE ****************************#

    def runDStatistic(self):
        try:
            self.statisticsCalculations.dAlignment = self.checkEntryPopulated(self.dAlignmentEntry, errorTitle='Missing Alignment', errorMessage='Please select and alignment.')
            self.statisticsCalculations.dWindowSize = self.checkEntryInRange(self.dWindowSizeEntry, min=0, inclusive=False, errorTitle='Invalid Window Size', errorMessage='Window size needs to be a positive integer.')
            self.statisticsCalculations.dWindowOffset = self.checkEntryInRange(self.dWindowOffsetEntry, min=0, inclusive=False, errorTitle='Invalid Window Offset', errorMessage='Window offset needs to be a positive integer.')
            self.statisticsCalculations.taxons = [self.dTaxonComboBox1.currentText(), self.dTaxonComboBox2.currentText(), self.dTaxonComboBox3.currentText(), self.dTaxonComboBox4.currentText()]

        except ValueError, (ErrorTitle, ErrorMessage, ErrorDescription):
            self.message(str(ErrorTitle), str(ErrorMessage), str(ErrorDescription))
            return

        self.statisticsCalculations.start()

    def displayDStatistic(self, dVal, dWindows):
        self.dVal = dVal
        self.dWindows = dWindows
        self.dStatisticWindow = dStatisticWindow.DStatisticWindow(self.dWindows)

        self.dStatisticValueLabel.setText(str(self.dVal))
        self.dStatisticLabel.setEnabled(True)
        self.dStatisticValueLabel.setEnabled(True)

    # **************************** MS PAGE ****************************#

    # LEO - added a comment, look at me!

    additionalFileCounter = 0
    additionalFileEntryNames = []

    def runMSCompare(self):
        try:
            self.msComparison.msToRax = False
            self.msComparison.msFiles = []
            self.msComparison.msTruth = self.checkEntryPopulated(self.msFileEntry, errorTitle='Missing MS Truth File', errorMessage='Please select an MS Truth file.')

            if self.checkboxCompareAgainstMS.isChecked():
                self.msComparison.msFiles.append(self.msSecondFileEntry.text())

                for i in range(len(self.additionalFileEntryNames)):
                    entry = self.findChild(QtGui.QLineEdit, self.additionalFileEntryNames[i])
                    self.msComparison.msFiles.append(self.checkEntryPopulated(entry, errorTitle='Blank Field', errorMessage='Field ' + str(i + 1) + ' is blank. Please select a file.'))

            if self.checkboxCompareAgainstRaxml.isChecked():
                self.msComparison.msToRax = True

                self.msComparison.raxmlDir = self.checkEntryPopulated(self.msRaxmlDirectoryEntry)
                self.msComparison.windowSize = int(self.checkEntryPopulated(self.msWindowSizeEntry))
                self.msComparison.windowOffset = int(self.checkEntryPopulated(self.msWindowOffsetEntry))

            self.msComparison.robinsonFouldsBarPlot = self.checkboxRobinsonFouldsBarPlot.isChecked()
            self.msComparison.percentMatchingSitesBarPlot = self.checkboxPercentMatchingSitesGraph.isChecked()
            self.msComparison.tmrcaLineGraph = self.checkboxTMRCAGraph.isChecked()

            if not (self.checkboxCompareAgainstRaxml.isChecked() or self.checkboxCompareAgainstMS.isChecked()):
                raise ValueError('Nothing to Compare Against', 'Please compare against a raxml directory and/or additional MS files.', 'n/a')

            if not (self.checkboxRobinsonFouldsBarPlot.isChecked() or self.checkboxPercentMatchingSitesGraph.isChecked() or self.checkboxTMRCAGraph.isChecked()):
                raise ValueError('No Plots Selected', 'Please select at least one plot.', 'n/a')

            self.msComparison.start()

        except ValueError, (ErrorTitle, ErrorMessage, ErrorDescription):
            self.message(ErrorTitle, ErrorMessage, ErrorDescription)
            return

    def plotMSCompare(self, unweightedData, percentMatchingSitesUnweighted, sitesToNewickMsMaps, msFiles, msTruthLabel):
        if self.msComparison.robinsonFouldsBarPlot:
            self.msRobinsonFouldsWindow = msRobinsonFouldsWindow.MSRobinsonFouldsWindow('Unweighted', unweightedData)

        if self.msComparison.percentMatchingSitesBarPlot:
            self.msPercentMatchingWindow = msPercentMatchingWindow.MSPercentMatchingWindow('Unweighted', percentMatchingSitesUnweighted, groupLabels1=msFiles)

        if self.msComparison.tmrcaLineGraph:
            self.msTMRCAWindow = msTMRCAWindow.MSTMRCAWindow(sitesToNewickMsMaps, [msTruthLabel] + msFiles)

    def addFileEntry(self, horizontalLayoutName, entryName, btnName, btn2Name):
        self.additionalFileCounter += 1
        self.additionalFileEntryNames.append(entryName + str(self.additionalFileCounter))

        # create horizontal layout
        HL = QtGui.QHBoxLayout()
        HL.setObjectName(horizontalLayoutName + str(self.additionalFileCounter))

        # create btn and add to horizontal layout
        btn2 = QtGui.QToolButton(self.msMSCompareGroupBox)
        btn2.setObjectName(btn2Name + str(self.additionalFileCounter))
        btn2.setText('-')
        btn2.setFixedHeight(21)
        btn2.setFixedWidth(23)
        HL.addWidget(btn2)

        # create text entry and add to horizontal layout
        entry = QtGui.QLineEdit(self.msMSCompareGroupBox)
        entry.setReadOnly(True)
        entry.setObjectName(entryName + str(self.additionalFileCounter))
        HL.addWidget(entry)

        # create btn and add to horizontal layout
        btn = QtGui.QToolButton(self.msMSCompareGroupBox)
        btn.setObjectName(btnName + str(self.additionalFileCounter))
        btn.setText('...')
        HL.addWidget(btn)

        self.resize(self.width(), self.height() + 30)
        self.msFileUploadMasterVL.addLayout(HL)

        btn.clicked.connect(lambda: self.getFileName(entry))
        btn2.clicked.connect(lambda: self.removeFileEntry(HL, entry, btn, btn2))

    def removeFileEntry(self, HL, entry, btn, btn2):
        HL.deleteLater()
        entry.deleteLater()
        btn.deleteLater()
        btn2.deleteLater()
        self.additionalFileEntryNames.remove(entry.objectName())
        self.resize(self.width(), self.height() - 30)

    # **************************** CONVERTER PAGE ****************************#

    def convertFile(self):
        try:
            self.fileToBeConverted = self.checkEntryPopulated(self.fileConverterEntry, errorTitle='No Input File Selected', errorMessage='Please choose an input file.')
            self.convertedFileDirectory = self.checkEntryPopulated(self.outputFileConverterEntry, errorTitle='No Output File Selected', errorMessage='Please choose an output file.')
        except ValueError, (ErrorTitle, ErrorMessage, ErrorDescription):
            self.message(ErrorTitle, ErrorMessage, ErrorDescription)
            return

        self.convertedFileName = self.convertedFileDirectory + '/convertedFile.' + self.outputFormatComboBox.currentText().lower() + '.txt'

        self.fileConverter.inputFileName = self.fileToBeConverted
        self.fileConverter.outputFileName = self.convertedFileName
        self.fileConverter.inputFormat = self.inputFormatComboBox.currentText().lower()
        self.fileConverter.outputFormat = self.outputFormatComboBox.currentText().lower()
        self.fileConverter.start()

    # **************************** RAXML PAGE ****************************#

    def generateSpeciesTree(self):
        try:
            # get values from gui -- ensure that no fields are blank
            self.raxmlOperations.inputFilename = self.checkEntryPopulated(self.inputFileEntry, errorTitle='Missing Alignment', errorMessage='Please select an alignment.')
            self.raxmlOperations.windowSize = self.checkEntryInRange(self.windowSizeEntry, min=0, inclusive=False, errorTitle='Invalid Window Size', errorMessage='Window size needs to be a positive integer.')
            self.raxmlOperations.windowOffset = self.checkEntryInRange(self.windowOffsetEntry, min=0, inclusive=False, errorTitle='Invalid Window Offset', errorMessage='Window offset needs to be a positive integer.')
            self.raxmlOperations.speciesTreeRooted = self.checkboxSpeciesTreeRooted.isChecked()
            self.raxmlOperations.speciesTreeOutGroup = self.speciesTreeComboBox.currentText()
            self.raxmlOperations.speciesTreeUseCustomRax = self.checkboxSpeciesTreeUseCustomRax.isChecked()

            # if using custom rax -- make sure that the user doesn't use the -s or -n flags
            self.raxmlOperations.speciesTreeCustomRaxmlCommand = ''
            if self.checkboxSpeciesTreeUseCustomRax.isChecked():
                self.raxmlOperations.speciesTreeCustomRaxmlCommand = self.checkEntryPopulated(self.speciesTreeRaxmlCommandEntry, errorTitle='No RAxML Command', errorMessage='Please enter a custom raxml command or uncheck the box.')
                if re.search('([\-][n])|([\-][s])', self.speciesTreeRaxmlCommandEntry.text()):
                    raise ValueError('Invalid RAxML Command', 'Please do not specify the -s or -n flags.', 'the -s and -n flags will be handled internally based on the alignment you input.')

        except ValueError, (ErrorTitle, ErrorMessage, ErrorDescription):
            self.message(str(ErrorTitle), str(ErrorMessage), str(ErrorDescription))
            return

        self.raxmlOperations.raxml_species_tree(self.raxmlOperations.inputFilename, rooted=self.raxmlOperations.speciesTreeRooted, outgroup=self.raxmlOperations.speciesTreeOutGroup, customRax=self.raxmlOperations.speciesTreeUseCustomRax, customRaxCommand=self.raxmlOperations.speciesTreeCustomRaxmlCommand)

    def updatedDisplayWindows(self, btnClicked=None):
        if btnClicked == None or btnClicked.isChecked():
            if self.runComplete == True:
                if self.raxmlInputErrorHandling():
                    # run commands that are shared by all functions
                    if self.getNumberChecked() > 0:
                        num = self.topTopologies
                        topologies_to_counts, unique_topologies_to_newicks = self.topologyPlotter.topology_counter(rooted=self.rooted, outgroup=self.outgroupComboBox.currentText())
                        self.numberOfUniqueTopologiesLabel.setText(str(len(topologies_to_counts)))
                        if num > len(topologies_to_counts):
                            num = len(topologies_to_counts)
                        list_of_top_counts, labels, sizes = self.topologyPlotter.top_freqs(num, topologies_to_counts)
                        top_topologies_to_counts = self.topologyPlotter.top_topologies(num, topologies_to_counts)
                        windows_to_top_topologies, top_topologies_list = self.topologyPlotter.windows_to_newick(top_topologies_to_counts, unique_topologies_to_newicks, rooted=self.rooted, outgroup=self.outgroupComboBox.currentText())  # all trees, scatter, circle, donut
                        topologies_to_colors, scatter_colors, ylist = self.topologyPlotter.topology_colors(windows_to_top_topologies, top_topologies_list)  # scatter, circle, (donut?)

                    # generate robinson foulds and pgtst graphs
                    if (btnClicked == None and self.checkboxRobinsonFoulds.isChecked()) or btnClicked == self.checkboxRobinsonFoulds:
                        if self.checkboxWeighted.isChecked():
                            windows_to_w_rf, windows_to_uw_rf = self.statisticsCalculations.calculate_windows_to_rf(self.speciesTree, self.checkboxWeighted.isChecked())
                            self.robinsonFouldsWindow = robinsonFouldsWindow.RobinsonFouldsWindow('Weighted Robinson Foulds Distance', windows_to_w_rf, 'Unweighted Robinson Foulds Distance', windows_to_uw_rf)
                        else:
                            windows_to_uw_rf = self.statisticsCalculations.calculate_windows_to_rf(self.speciesTree, self.checkboxWeighted.isChecked())
                            self.robinsonFouldsWindow = robinsonFouldsWindow.RobinsonFouldsWindow('Weighted Robinson Foulds Distance', windows_to_uw_rf)

                    if (btnClicked == None and self.checkboxPGTST.isChecked()) or btnClicked == self.checkboxPGTST:
                        windowsToPGTST = self.statisticsCalculations.calculate_windows_to_p_gtst(self.speciesTree)
                        self.pgtstWindow = pgtstWindow.PGTSTWindow(windowsToPGTST, "p(gt|st)", xLabel="Windows", yLabel="Probability")

                    # generate donut plot
                    if (btnClicked == None and self.checkboxDonutPlot.isChecked()) or btnClicked == self.checkboxDonutPlot:
                        donut_colors = self.topologyPlotter.donut_colors(top_topologies_to_counts, topologies_to_colors)  # donut
                        self.donutPlotWindow = donutPlotWindow.DonutPlotWindow('Frequency of Top Topologies', labels, sizes, donut_colors)

                    # generate scatter plot
                    if (btnClicked == None and self.checkboxScatterPlot.isChecked()) or btnClicked == self.checkboxScatterPlot:
                        self.scatterPlotWindow = scatterPlotWindow.ScatterPlotWindow('Windows to Top Topologies', windows_to_top_topologies, scatter_colors, ylist)

                    # generate circle graph
                    if (btnClicked == None and self.checkboxCircleGraph.isChecked()) or btnClicked == self.checkboxCircleGraph:
                        sites_to_informative, windows_to_informative_count, windows_to_informative_pct, pct_informative = self.informativeSites.calculate_informativeness('windows', self.raxmlOperations.windowOffset)
                        self.circleGraphWindow = circleGraphWindow.CircleGraphWindow(self.raxmlOperations.inputFilename, windows_to_top_topologies, topologies_to_colors, self.raxmlOperations.windowSize, self.raxmlOperations.windowOffset, sites_to_informative)

                    # generate informative sites heatmap graph
                    if (btnClicked == None and self.checkboxHeatMap.isChecked()) or btnClicked == self.checkboxHeatMap:
                        sites_to_informative, windows_to_informative_count, windows_to_informative_pct, pct_informative = self.informativeSites.calculate_informativeness('windows', self.raxmlOperations.windowOffset)
                        self.heatMapWindow = heatMapWindow.HeatMapWindow('Heat Map', sites_to_informative)

                    # generate windows to informative sites line graph
                    if (btnClicked == None and self.checkboxWindowsToInfSites.isChecked()) or btnClicked == self.checkboxWindowsToInfSites:
                        sites_to_informative, windows_to_informative_count, windows_to_informative_pct, pct_informative = self.informativeSites.calculate_informativeness('windows', self.raxmlOperations.windowOffset)
                        self.windowsToInfSitesWindow = windowsToInfSitesWindow.WindowsToInfSitesWindow('Windows to Informative Sites', windows_to_informative_pct)

                    # generate bootstrap graph
                    if (btnClicked == None and self.checkboxBootstrap.isChecked()) or btnClicked == self.checkboxBootstrap:
                        internal_nodes_i, internal_nodes_f = self.bootstrapContraction.internal_nodes_after_contraction(self.confidenceLevel)
                        self.bootstrapContractionWindow = bootstrapContractionWindow.BootstrapContractionWindow(internal_nodes_i, internal_nodes_f, self.confidenceLevel, xLabel="Window Indices", yLabel="Number of Internal Nodes")

                    # generate all trees graph
                    if (btnClicked == None and self.checkboxAllTrees.isChecked()) or btnClicked == self.checkboxAllTrees:
                        self.allTreesWindow = allTreesWindow.AllTreesWindow('', topologies_to_colors, rooted=self.checkboxRooted.isChecked(), outGroup=self.outgroupComboBox.currentText())

    def raxmlInputErrorHandling(self):
        """
            returns true if all tests pass otherwise false
        """
        try:
            # input alignment for raxml
            self.raxmlOperations.inputFilename = self.checkEntryPopulated(self.inputFileEntry, errorTitle='Missing Alignment', errorMessage='Please select an alignment.')
            self.raxmlOperations.windowSize = self.checkEntryInRange(self.windowSizeEntry, min=0, inclusive=False, errorTitle='Invalid Window Size', errorMessage='Window size needs to be a positive integer.')
            self.raxmlOperations.windowOffset = self.checkEntryInRange(self.windowOffsetEntry, min=0, inclusive=False, errorTitle='Invalid Window Offset', errorMessage='Window offset needs to be a positive integer.')
            self.raxmlOperations.outGroup = self.outgroupComboBox.currentText()
            self.raxmlOperations.model = self.modelComboBox.currentText()
            self.raxmlOperations.isCustomRaxmlCommand = self.checkBoxCustomRaxml.isChecked()
            self.raxmlOperations.bootstrap = self.checkboxBootstrap.isChecked()
            self.raxmlOperations.rooted = self.checkboxRooted.isChecked()
            self.rooted = self.checkboxRooted.isChecked()

            # if user is generating Top Topologies or scatter plot or donut plor or circle graph run error handling on top topologies entry
            if self.checkboxAllTrees.isChecked() or self.checkboxScatterPlot.isChecked() or self.checkboxDonutPlot.isChecked() or self.checkboxCircleGraph.isChecked():
                self.checkEntryPopulated(self.numberOfTopTopologiesEntry, errorTitle='Number of Top Topologies Field is Blank', errorMessage='Please enter a number of top topologies.')
                self.topTopologies = self.checkEntryInRange(self.numberOfTopTopologiesEntry, min=0, max=16, inclusive=False, errorTitle='Invalid Number of Top Topologies', errorMessage='Please enter an integer between 0 and 15.')

            # bootstrap error handling
            self.raxmlOperations.numBootstraps = 0
            if self.checkboxBootstrap.isChecked():
                self.confidenceLevel = self.checkEntryInRange(self.confidenceLevelEntry, min=0, max=100, errorTitle='Invalid Confidence Level', errorMessage='Please enter an integer between 0 and 100.')
                self.raxmlOperations.numBootstraps = self.checkEntryInRange(self.numberOfBootstrapsEntry, min=2, errorTitle='Invalid Number of Bootstraps', errorMessage='Please enter an integer greater than 1.')

            # if using custom rax -- make sure that the user doesn't use the -s or -n flags
            if self.checkBoxCustomRaxml.isChecked():
                self.raxmlOperations.customRaxmlCommand = self.checkEntryPopulated(self.customRaxmlCommandEntry, errorTitle='No RAxML Command', errorMessage='Please enter a custom raxml command or uncheck the box.')
                if re.search('([\-][n])|([\-][s])', self.customRaxmlCommandEntry.text()):
                    raise ValueError, ('Invalid RAxML Command', 'Please do not specify the -s or -n flags.', 'the -s and -n flags will be handled internally based on the alignment you input.')

            # species tree error handling
            if self.speciesTreeEntry.text() != "" and self.newickFileEntry.text() != "":
                raise ValueError, ('Multiple Species Trees', 'You have both selected a species tree file and entered a species tree. Please only do one.', 'Both the "Species Tree File and "Enter Species Tree" fields are populated. Please only use one.')

            # if the user selects either statistic plot -- open the inputted newick and read it into memory as a string on a single line
            if self.checkboxRobinsonFoulds.isChecked() or self.checkboxPGTST.isChecked():
                if self.newickFileEntry.text() != "":
                    self.newickFileName = self.checkEntryPopulated(self.newickFileEntry, errorTitle='Missing Species Tree', errorMessage='Please select a species tree.', errorDescription='Please select a species tree.')
                    with open(self.newickFileEntry.text(), 'r') as f:
                        self.speciesTree = f.read().replace('\n', '')
                else:
                    self.speciesTree = self.checkEntryPopulated(self.speciesTreeEntry, errorTitle='Missing Species Tree', errorMessage='Please select a species tree.', errorDescription='Please select a species tree.')


        except ValueError, (ErrorTitle, ErrorMessage, ErrorDescription):
            self.message(str(ErrorTitle), str(ErrorMessage), str(ErrorDescription))
            return False

        return True

    def runRAxML(self):
        # if all error handling passes run RAxML
        if self.raxmlInputErrorHandling():
            if self.runComplete:
                self.msg = QtGui.QMessageBox()
                self.msg.setText("Rerun RAxML?")
                self.msg.setInformativeText('Are you sure you want to rerun RAxML?')
                self.msg.setStandardButtons(QtGui.QMessageBox.Yes | QtGui.QMessageBox.Cancel)

                # set icon
                pixmap = QtGui.QPixmap('imgs/warning.png')
                self.msg.setIconPixmap(pixmap)

                # execute window
                returnVal = self.msg.exec_()

                # if the user selected the 'ok' button
                if returnVal == QtGui.QMessageBox.Yes:
                    # start raxml operations thread
                    self.raxmlOperations.start()
            # if raxml hasn't been run before just run it
            else:
                # start raxml operations thread
                self.raxmlOperations.start()

    def raxmlComplete(self):
        topologies_to_counts, unique_topologies_to_newicks = self.topologyPlotter.topology_counter(rooted=self.rooted, outgroup=self.outgroupComboBox.currentText())
        self.numberOfUniqueTopologiesLabel.setText(str(len(topologies_to_counts)))
        self.progressBar.setValue(100)
        self.runComplete = True

    # **************************** ABSTRACT ****************************#

    def message(self, title, description, extraInfo, type='Err'):
        """
            creates and displays and window displaying the message
        """

        # create object
        errMessage = QtGui.QMessageBox()

        # set text
        errMessage.setText(title)
        errMessage.setInformativeText(description)
        errMessage.setDetailedText(extraInfo)

        # default pixmap for error
        pixmap = QtGui.QPixmap('imgs/warning.png')

        # set icon
        errMessage.setIconPixmap(pixmap)

        # execute window
        errMessage.exec_()

    def checkEntryPopulated(self, entry, errorTitle='Field Not Populated', errorMessage='Please populate field.', errorDescription=None):
        """
            checks if given entry is empty or not.
                (i) if entry is populated returns text
                (ii) otherwise raises value error
        """

        # if user does not provide an error description generate one automatically
        if not errorDescription:
            errorDescription = 'relevant entry name: ' + str(entry.objectName())

        text = str(entry.text())

        if text == '':
            raise ValueError(errorTitle, errorMessage, errorDescription)

        return text

    def checkEntryInRange(self, entry, min=(-1.0 * float('inf')), max=float('inf'), inclusive=True, errorTitle='Entry Out Of Range', errorMessage='', errorDescription=None):
        """
            checks if value of given entry is in range.
                i. if entry is in given range return it
                ii. otherwise raises value error
        """

        # if user does not provide an error description generate one automatically
        if not errorDescription:
            errorDescription = 'relevant entry name: ' + str(entry.objectName())

        # check to make sure the entry is populated
        if entry.text() != '':
            val = float(int(float(entry.text())))
        else:
            raise ValueError, (errorTitle, errorMessage, errorDescription)

        # check to make sure value is in range
        if inclusive:
            if val < min or val > max:
                raise ValueError, (errorTitle, errorMessage, errorDescription)
        else:
            if val <= min or val >= max:
                raise ValueError, (errorTitle, errorMessage, errorDescription)

        return int(val)

    def updateTaxonComboBoxes(self, comboBoxes, textEntry, require4Taxons=False):
        """
            input:
                i. comboBoxes - a list of comboBox widgets (drop down menus)
                ii. textEntry - a text entry widget
                iii. errHandling=False - a boolean indicating whether or not to require that there be exactly four taxons in the file in the text entry

            gets a list of taxons from the file in textEntry and sets the items in a list of combo boxes to that list of taxons.
        """
        try:
            if textEntry.text() == "":
                return

            # get list of taxon names from file
            taxonNames = list(self.raxmlOperations.taxon_names_getter(textEntry.text()))

            if require4Taxons:
                # if there are not exactly 4 taxons
                if len(taxonNames) != 4:
                    self.message('Invalid File.', 'Need exactly 4 taxons.', textEntry.text())
                    return

            # clear each combo box
            for comboBox in comboBoxes:
                comboBox.clear()

            # add the list of taxons to each combobox
            for taxon in taxonNames:
                for comboBox in comboBoxes:
                    comboBox.addItem(taxon)

            for i in range(len(comboBoxes)):
                comboBoxes[i].setCurrentIndex(i)

        except:
            self.message('Invalid File', 'Invalid alignment file. Please choose another.', 'Make sure your file has 4 sequences and is in the phylip-relaxed format.', type='Err')
            return

    def getNumberChecked(self):
        """
            returns the number of checkboxes that are checked
        """
        return (self.checkboxScatterPlot.checkState() + self.checkboxCircleGraph.checkState() + self.checkboxDonutPlot.checkState() + self.checkboxAllTrees.checkState()) / 2

    def toggleEnabled(self, guiElement):
        """
           toggles whether or not guiElement is enabled
        """
        enabled = guiElement.isEnabled()
        guiElement.setEnabled(not enabled)

    def setWindow(self, window):
        self.stackedWidget.setCurrentIndex(self.windows[window])
        self.resize(self.windowSizes[window]['x'], self.windowSizes[window]['y'])
        self.move(self.windowLocations[window]['x'], self.windowLocations[window]['y'])

    def ensureSingleModeSelected(self, mode_selected, window):
        for mode in self.menuMode.actions():
            if mode != mode_selected:
                mode.setChecked(False)

        mode_selected.setChecked(True)
        self.setWindow(window)

    def saveFileAs(self, textEntry):
        """
            i. open a dialog to get in which user enters a file name to save
            ii. sets the text of given text entry to match file user selected
        """
        textEntry.setText(QtGui.QFileDialog.getSaveFileName(self, 'Export'))

    def getFileName(self, textEntry):
        """
            i. open a dialog to get in which user selects a file
            ii. sets the text of given text entry to match file user selected
        """
        textEntry.setText(QtGui.QFileDialog.getOpenFileName())
        textEntry.emit(QtCore.SIGNAL('FILE_SELECTED'))

    def openDirectory(self, textEntry):
        """
            i. open a dialog in which user selects a directory
            ii. sets the text of given text entry to match the directory the user selected
        """
        textEntry.setText(QtGui.QFileDialog.getExistingDirectory())
        textEntry.emit(QtCore.SIGNAL("DIRECTORY_SELECTED"))

    def openWindow(self, window, type='std'):
        window.show()
        if type == 'std':
            window.plot()
        elif type == 'tabs':
            window.displayImages()

    def openURL(self, url):
        webbrowser.open(url, new=0, autoraise=True)

    def setValidator(self, entry, validator):
        if validator == 'Double':
            entry.setValidator(QtGui.QDoubleValidator(entry))
        elif validator == 'Int':
            entry.setValidator(QtGui.QIntValidator(entry))

    # def resizeEvent(self, event):
    #     print self.size()

    # def moveEvent(self, QMoveEvent):
    #     print self.pos()


if __name__ == '__main__':  # if we're running file directly and not importing it
    app = QtGui.QApplication(sys.argv)  # A new instance of QApplication

    # initialize main input window
    form = PhyloVisApp()  # We set the form to be our PhyloVisApp (design)
    form.show()  # Show the form
    form.move(600, 300)

    sys.exit(app.exec_())  # and execute the app
