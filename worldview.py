from PyQt4 import QtCore, QtGui, QtOpenGL
import pyqtgraph.opengl as gl
import numpy as np

app = QtGui.QApplication([])    #Starts up the QtGui; makes gl object that needs to be passed to graph routines...

class WorldView(gl.GLViewWidget):
    def __init__(self, parent=None):
        #global ShareWidget

        #if ShareWidget is None:
        #    # create a dummy widget to allow sharing objects
        #    # (textures, shaders, etc) between views
        ShareWidget = QtOpenGL.QGLWidget()

        QtOpenGL.QGLWidget.__init__(self, parent, ShareWidget)

        self.setFocusPolicy(QtCore.Qt.ClickFocus)

        self.opts = {
            # will always appear at the center of the widget
            #'center': Vector(0, 0, 0),
            'center': [0,0,0],

            # distance of camera from center
            'distance': 10.0,
            # horizontal field of view in degrees
            'fov':  60,
            # camera's angle of elevation in degrees
            'elevation':  30,
            # camera's azimuthal angle in degrees
            # (rotation around z-axis 0 points along x-axis)
            'azimuth': 45,
            # glViewport params; None == whole widget
            'viewport': None,
        }
        self.setBackgroundColor('k')
        self.items = []
        self.noRepeatKeys = [QtCore.Qt.Key_Right, QtCore.Qt.Key_Left,
                             QtCore.Qt.Key_Up, QtCore.Qt.Key_Down,
                             QtCore.Qt.Key_PageUp, QtCore.Qt.Key_PageDown]
        self.keysPressed = {}
        self.keyTimer = QtCore.QTimer()
        self.keyTimer.timeout.connect(self.evalKeyState)

        self.makeCurrent()

        self.ray = QtGui.QVector3D(0, 0, 0)
        self.select = False

    def mousePressEvent(self, event):
        print ("Pressed button", event.button(), "at", event.pos())

        self.mousePos = event.pos()
        if event.button() == 2:
            self.select = True
        else:
            self.select = False
        print (self.itemsAt((self.mousePos.x(), self.mousePos.y(), 3, 3)))


    def mousePos(self, ev):
        """ Store the position of the mouse press for later use.

        """
        self.mousePos = ev.pos()

b = []
for k in range(250000):
    b.append(np.random.rand(3))
b=np.array(b)

w = WorldView()
sp1 = gl.GLScatterPlotItem(pos=b, size=1)
w.addItem(sp1)

__WIDTH = 512
__HEIGHT = 424

m = w.projectionMatrix() * w.viewMatrix()

projected_array = np.zeros((__WIDTH * __HEIGHT, 2))
view_w = w.width()
view_h = w.height()
w.mousePos(w)
    
mouse_x = w.mousePos.x()
mouse_y = w.mousePos.y()

# b array contains the raw coordinates of all the points on screen      
step=1  
for i in xrange(0, __WIDTH, step):
    for j in xrange(0, __HEIGHT, step):
        pt = m.map(QtGui.QVector3D(b[j*__WIDTH+i, 0], 
                                   b[j*__WIDTH+i, 1], 
                                   b[j*__WIDTH+i, 2]))
        
        # origin range [-1, 1]
        projected_array[j*__WIDTH+i, 0] = (pt.x() + 1)/2
        projected_array[j*__WIDTH+i, 1] = (- pt.y() + 1)/2


projected_array[:, 0] = (projected_array[:, 0] - (mouse_x/view_w))
projected_array[:, 1] = (projected_array[:, 1] - (mouse_y/view_h))
distance_array = np.power(np.power(projected_array[:, 0], 2) +
                          np.power(projected_array[:, 1], 2), 0.5)

min_index = np.nanargmin(distance_array)
print min_index
# mark the selected point on screen with a big sphere
#sp2.setData(pos=b[min_index, :], size=100)


### Start Qt event loop unless running in interactive mode.
#QtGui.QApplication.instance().exec_()
#app.closeAllWindows()
