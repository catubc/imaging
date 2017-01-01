# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'TestWindow.ui'
#
# Created: Fri Jan 17 11:13:29 2014
#      by: PyQt4 UI code generator 4.9.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(911, 768)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.openglwidget = gloverride(self.centralwidget)
        self.openglwidget.setGeometry(QtCore.QRect(389, 16, 270, 309))
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(100)
        sizePolicy.setVerticalStretch(100)
        sizePolicy.setHeightForWidth(self.openglwidget.sizePolicy().hasHeightForWidth())
        self.openglwidget.setSizePolicy(sizePolicy)
        self.openglwidget.setStyleSheet(_fromUtf8("border:none;"))
        self.openglwidget.setObjectName(_fromUtf8("openglwidget"))
        self.topgraphicsView = PlotWidget(self.centralwidget)
        self.topgraphicsView.setGeometry(QtCore.QRect(98, 138, 256, 227))
        self.topgraphicsView.setStyleSheet(_fromUtf8("border:none;"))
        self.topgraphicsView.setObjectName(_fromUtf8("topgraphicsView"))
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 911, 18))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "MainWindow", None, QtGui.QApplication.UnicodeUTF8))

from GLsubclass import gloverride
from pyqtgraph import PlotWidget
