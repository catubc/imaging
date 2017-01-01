import pyqtgraph as pg
import pyqtgraph.opengl

def addSphere(v, r, c):
    ms = pyqtgraph.opengl.MeshData.sphere(rows=10, cols=30, radius=r)
    gs = pyqtgraph.opengl.GLMeshItem(meshdata=ms, smooth=True, drawFaces=True, color=c, drawEdges=False, shader="shaded", glOptions=('opaque' if c[-1] == 1. else 'translucent'))
    gs.translate(v[0], v[1], v[2])
    view.addItem(gs)

class MyView(pg.opengl.GLViewWidget):
    def mousePressEvent(self, ev):
        print("Pressed button", ev.button(), "at", (ev.pos().x(), ev.pos().y()), " items: ", self.itemsAt((ev.pos().x(), ev.pos().y(), 10, 10)))

        # If you do not accept the event, then no move/release
        # events will be received.
        ev.accept()


pg.mkQApp()
view = MyView()
addSphere((0,0,0),1,(1,1,1,1))
view.show()

