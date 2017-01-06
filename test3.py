import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



data = np.random.random((100,3))
print data.shape


from scipy.spatial import ConvexHull
hull = ConvexHull(data)
print "...volume: ", hull.volume

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for k in range(0,len(data),1): 
    ax.scatter(data[k][0], data[k][1], data[k][2], marker='o')

for simplex in hull.simplices:
    plt.plot(data[simplex, 0], data[simplex, 1], data[simplex, 2] , 'k-')



ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()
            
            
            
