import os, sys, matplotlib, matplotlib.pyplot
import numpy as np
from pyqtgraph.Qt import QtCore, QtGui
from pyqtgraph.widgets.GraphicsLayoutWidget import GraphicsLayoutWidget
import pyqtgraph as pg
import pyqtgraph.functions as fn
import matplotlib.pyplot as plt

N = 256
ARR = np.random.random((N,N))*255
norm = plt.Normalize()
ARR_OFF = ARR #plt.cm.jet(norm(ARR))
# Change ARR_OFF to ARR to see my problem

class MainWindow(QtGui.QMainWindow):

    def __init__(self, parent=None):

        QtGui.QMainWindow.__init__(self, parent)
        self.setupUserInterface()
        self.setupSignals()

    def setupUserInterface(self):
        """ Initialise the User Interface """
        # Left frame
        leftFrame = QtGui.QFrame()
        leftFrameLayout = QtGui.QHBoxLayout()
        leftFrame.setLayout(leftFrameLayout)
        leftFrame.setLineWidth(0)
        leftFrame.setFrameStyle(QtGui.QFrame.Panel)
        leftFrameLayout.setContentsMargins(0,0,5,0)

        # Left frame contents
        self.viewMain = GraphicsLayoutWidget()  # A GraphicsLayout within a GraphicsView
        leftFrameLayout.addWidget(self.viewMain)
        self.viewMain.setMinimumSize(200,200)
        self.vb = MultiRoiViewBox(lockAspect=True,enableMenu=True)
        self.viewMain.addItem(self.vb)
        self.vb.enableAutoRange()

        # Right frame
        self.sidePanel = SidePanel(self)

        # UI window (containing left and right frames)
        UIwindow         = QtGui.QWidget(self)
        UIwindowLayout   = QtGui.QHBoxLayout()
        UIwindowSplitter = QtGui.QSplitter(QtCore.Qt.Horizontal)
        UIwindowLayout.addWidget(UIwindowSplitter)
        UIwindow.setLayout(UIwindowLayout)
        self.setCentralWidget(UIwindow)
        UIwindowSplitter.addWidget(leftFrame)
        UIwindowSplitter.addWidget(self.sidePanel)

        self.setMinimumSize(600,500)
        self.resize(self.minimumSize())

    def setupSignals(self):
        """ Setup signals """
        self.sidePanel.buttImageAdd.clicked.connect(self.showImage)

    def showImage(self,imageFilename):
        """ Shows image in main view """
        self.vb.showImage(ARR)

class ViewMode():
    def __init__(self,id,cmap):
        self.id   = id
        self.cmap = cmap
        self.getLookupTable()
    def getLookupTable(self):
        lut = [ [ int(255*val) for val in self.cmap(i)[:3] ] for i in xrange(256) ]
        lut = np.array(lut,dtype=np.ubyte)
        self.lut = lut

class MultiRoiViewBox(pg.ViewBox):

    def __init__(self,parent=None,border=None,lockAspect=False,enableMouse=True,invertY=False,enableMenu=True,name=None):
        pg.ViewBox.__init__(self,parent,border,lockAspect,enableMouse,invertY,enableMenu,name)
        self.img      = None
        self.NORMAL   = ViewMode(0,matplotlib.cm.gray)
        self.DEXA     = ViewMode(1,matplotlib.cm.jet)
        self.viewMode = self.NORMAL

    def showImage(self,arr):
        if arr==None:
            self.img = None
            return
        if self.img==None:
            self.img = pg.ImageItem(arr,autoRange=False,autoLevels=False)
            self.addItem(self.img)
        self.img.setImage(arr,autoLevels=False)
        self.updateView()

    def updateView(self):
        self.background.setBrush(fn.mkBrush(self.viewMode.lut[0]))
        self.background.show()
        if    self.img==None: return
        else: self.img.setLookupTable(self.viewMode.lut)


from pyqtgraph.Qt import QtCore,QtGui

class SidePanel(QtGui.QWidget):

    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self,parent)
        self.setMinimumWidth(250)
        self.buttMinimumSize = QtCore.QSize(36,36)
        self.setupImageToolbox()
        sidePanelLayout = QtGui.QVBoxLayout()
        sidePanelLayout.addWidget(self.imageToolbox)
        sidePanelLayout.setContentsMargins(0,0,0,0)
        self.setLayout(sidePanelLayout)

    def setupImageToolbox(self):
        # Image buttons
        self.buttImageAdd  = QtGui.QPushButton()
        imageButtons       = [self.buttImageAdd]
        for i in xrange(len(imageButtons)):
            image = imageButtons[i]
            image.setMinimumSize(self.buttMinimumSize)

        self.imageFileTools  = QtGui.QFrame()
        imageFileToolsLayout = QtGui.QHBoxLayout()
        self.imageFileTools.setLayout(imageFileToolsLayout)
        self.imageFileTools.setLineWidth(1)
        self.imageFileTools.setFrameStyle(QtGui.QFrame.StyledPanel)
        imageFileToolsLayout.addWidget(self.buttImageAdd)

        # Image Toolbox (containing imageFileList + imageFileList buttons)
        self.imageToolbox = QtGui.QFrame()
        self.imageToolbox.setLineWidth(2)
        self.imageToolbox.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
        imageToolboxLayout = QtGui.QVBoxLayout()
        self.imageToolbox.setLayout(imageToolboxLayout)
        imageToolboxLayout.addWidget(self.imageFileTools)


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
