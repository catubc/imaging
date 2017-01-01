import pyqtgraph as pg
import pyqtgraph.opengl
from doubleGUI import Ui_MainWindow
from PyQt4 import QtCore, QtGui
import sys


class mainwindow(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None):
        super(mainwindow, self).__init__(parent=parent)

        self.setupUi(self)
        
def main():
    app = QtGui.QApplication(sys.argv)
    
    mywindow = mainwindow()
    mywindow.show()  
    sys.exit(app.exec_())
    
main()
