from pyqtgraph.opengl import GLViewWidget

class gloverride(GLViewWidget):
    def __init__(self, parent=None):
        print "bla."
    
        #self.mousePressEvent(ev)
        
    def mousePressEvent(self, ev):
        self.mousePos = ev.pos()
        print "Pressed button HERE", ev.button(), "at", ev.pos()


    def mouseReleaseEvent(self, ev):
        region = (ev.pos().x()-5, ev.pos().y()-5, 10, 10)
        # y inversion
        region = (region[0], self.height()-(region[1]+region[3]), region[2], region[3])

        # which items are here?
        print(self.itemsAt(region))

        # draw the picking region
        glViewport(*self.getViewport())
        glClear( GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT )
        self.paintGL(region=region)
        self.swapBuffers()
