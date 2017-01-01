import pyqtgraph as pg
import numpy as np
app = pg.mkQApp()

plt = pg.plot()

lines = 100
points = 20000
x = np.empty((lines, points))
x[:] = np.arange(points)
y = np.random.normal(size=(lines, points))
for k in range(len(y)):
    y[k]+=k*10
connect = np.ones((lines, points), dtype=np.ubyte)
connect[:,-1] = 0  #  disconnect segment between lines

path = pg.arrayToQPath(x.reshape(lines*points), y.reshape(lines*points), connect.reshape(lines*points))

item = pg.QtGui.QGraphicsPathItem(path)

item.setPen(pg.mkPen('w'))

plt.addItem(item)


app.exec_()
