import filter
import os
import glob
import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt
import struct
import string
import scipy
import time
import matplotlib.animation as animation
from matplotlib.path import Path

from threading import *
import time

import sklearn  
from sklearn import manifold
from tsne import bh_sne
import scipy.ndimage as ndimage

from pyqtgraph.Qt import QtCore, QtGui, QtOpenGL

import pyqtgraph as pg
import pyqtgraph.opengl as gl
import matplotlib.gridspec as gridspec

import tifffile as tiff
from matplotlib.colors import LinearSegmentedColormap


cdict1 = {'red':   ((0.0, 0.0, 0.0),
                   (0.5, 0.0, 0.1),
                   (1.0, 1.0, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.0, 1.0),
                   (0.5, 0.1, 0.0),
                   (1.0, 0.0, 0.0))
        }
        
blue_red1 = LinearSegmentedColormap('BlueRed1', cdict1)

def PCA(X, n_components):
    from sklearn import decomposition

    pca = decomposition.PCA(n_components)
    pca.fit(X)
    X=pca.transform(X)

    coords = []
    for i in range(len(X)):
         coords.append([X[i][0], X[i][1], X[i][2]])
    
    return X, np.array(coords).T #THIS IS REDUNDANT... REDUCE IT

def dim_reduction(matrix_in, method, mouse):
    
    methods = ['MDS - SMACOF', 't-SNE', 'PCA', 't-SNE Barnes Hut','not implemented','not implemented']
    print "Computing dim reduction, size of array: ", np.array(matrix_in).shape
        
    if method==0:
        #MDS Method - SMACOF implementation Nelle Varoquaux
        if os.path.exists(mouse+'_'+str(len(matrix_in))+'_MDS.npy')==False:
            print "... MDS-SMACOF..."
            print "... pairwise dist ..."
            dists = sklearn.metrics.pairwise.pairwise_distances(matrix_in.astype(np.float64))
            adist = np.array(dists)
            amax = np.amax(adist)
            adist /= amax
            
            print "... computing MDS ..."
            mds_clf = manifold.MDS(n_components=3, metric=True, n_jobs=-1, dissimilarity="precomputed", random_state=6)
            results = mds_clf.fit(adist)
            Y = results.embedding_ 

            np.save(mouse+'_'+str(len(matrix_in))+'_MDS', Y)
        else: print "...already computed...."
        #    Y = np.load(mouse+'_'+str(len(matrix_in))+'_MDS.npy')
                
    elif method==1:
        ##t-Distributed Stochastic Neighbor Embedding; Laurens van der Maaten
        if os.path.exists(mouse+str(len(matrix_in))+'_tSNE.npy')==False:
            print "... tSNE ..."
            print "... pairwise dist ..."
            
            dists = sklearn.metrics.pairwise.pairwise_distances(matrix_in)
            
            adist = np.array(dists)
            amax = np.amax(adist)
            adist /= amax
            
            print "... computing tSNE ..."
            model = manifold.TSNE(n_components=3, init='pca', random_state=0)
            Y = model.fit_transform(adist)
            #Y = model.fit(adist)
        
            np.save(mouse+'_'+str(len(matrix_in))+'_tSNE', Y)
        #else: print "...already computed...."
        #    Y = np.load(home_dir+mouse+'_tSNE.npy')

    elif method==2:

        if os.path.exists(mouse+'_'+str(len(matrix_in))+'_PCA.npy')==False:
            print "...computing PCA..."
            Y, X = PCA(matrix_in, 3)

            np.save(mouse+'_'+str(len(matrix_in))+'_PCA', Y)
        else: print "...already computed...."
        #    Y = np.load(mouse+'_'+str(len(matrix_in))+'_PCA.npy')
            
                
    elif method==3:

        if os.path.exists(home_dir+mouse+str(len(matrix_in))+'_tSNE_barnes_hut.npy')==False:
            print "... computing Barnes-Hut tSNE..."
            Y = bh_sne(np.array(matrix_in))
            
            np.save(home_dir+mouse+'_tSNE_barnes_hut', Y)
        else: print "...already computed...."
        #    Y = np.load(home_dir+mouse+'_tSNE_barnes_hut.npy')
            
    
    elif method==4: 
        print "NOT IMPLEMENTED"
        
        sammon(matrix_in)    

class MyGLViewWidget(gl.GLViewWidget):
    """ Override GLViewWidget with enhanced behavior and Atom integration.
    
    """
    #: Fired in update() method to synchronize listeners.
    sigUpdate = QtCore.pyqtSignal()
    
    def mousePressEvent(self, ev):
        """ Store the position of the mouse press for later use.
        
        """
        super(MyGLViewWidget, self).mousePressEvent(ev)
        self._downpos = self.mousePos
            
    def mouseReleaseEvent(self, ev):
        """ Allow for single click to move and right click for context menu.
        
        Also emits a sigUpdate to refresh listeners.
        """
        super(MyGLViewWidget, self).mouseReleaseEvent(ev)
        if self._downpos == ev.pos():
            if ev.button() == 2:
                print 'show context menu'
            elif ev.button() == 1:
                x = ev.pos().x() - self.width() / 2
                y = ev.pos().y() - self.height() / 2
                self.pan(-x, -y, 0, relative=True)
                print self.opts['center']
        self._prev_zoom_pos = None
        self._prev_pan_pos = None
        self.sigUpdate.emit()
                
    def mouseMoveEvent(self, ev):
        """ Allow Shift to Move and Ctrl to Pan.
        
        """
        shift = ev.modifiers() & QtCore.Qt.ShiftModifier
        ctrl = ev.modifiers() & QtCore.Qt.ControlModifier
        if shift:
            y = ev.pos().y()
            if not hasattr(self, '_prev_zoom_pos') or not self._prev_zoom_pos:
                self._prev_zoom_pos = y
                return
            dy = y - self._prev_zoom_pos
            def delta():
                return -dy * 5
            ev.delta = delta
            self._prev_zoom_pos = y
            self.wheelEvent(ev)
        elif ctrl:
            pos = ev.pos().x(), ev.pos().y()
            if not hasattr(self, '_prev_pan_pos') or not self._prev_pan_pos:
                self._prev_pan_pos = pos
                return
            dx = pos[0] - self._prev_pan_pos[0]
            dy = pos[1] - self._prev_pan_pos[1]
            self.pan(dx, dy, 0, relative=True)
            self._prev_pan_pos = pos
        else:
            super(MyGLViewWidget, self).mouseMoveEvent(ev)


        
class WorldView(gl.GLViewWidget):
    #: Fired in update() method to synchronize listeners.
    sigUpdate = QtCore.pyqtSignal()
    
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


class MyGLViewWidget(gl.GLViewWidget):
    """ Override GLViewWidget with enhanced behavior and Atom integration.
    
    """
    #: Fired in update() method to synchronize listeners.
    sigUpdate = QtCore.pyqtSignal()
    
    def mousePressEvent(self, ev):
        global pts
        b = pts
        """ Store the position of the mouse press for later use.
        
        """
        #super(MyGLViewWidget, self).mousePressEvent(ev)
        #self._downpos = self.mousePos


        #__WIDTH = 512
        #__HEIGHT = 424
        
        #m = self.projectionMatrix() * self.viewMatrix()
        #projected_array = np.zeros((__WIDTH * __HEIGHT, 2))
        #view_w = self.width()
        #view_h = self.height()
        #mouse_x = self.mousePos.x()
        #mouse_y = self.mousePos.y()
        
        ## b array contains the raw coordinates of all the points on screen  
        #step=1      
        #for i in xrange(0, __WIDTH, step):
            #for j in xrange(0, __HEIGHT, step):
                #pt = m.map(QtGui.QVector3D(b[j*__WIDTH+i, 0],
                                           #b[j*__WIDTH+i, 1],
                                           #b[j*__WIDTH+i, 2]))
                ## origin range [-1, 1]
                #projected_array[j*__WIDTH+i, 0] = (pt.x() + 1)/2
                #projected_array[j*__WIDTH+i, 1] = (- pt.y() + 1)/2


        #projected_array[:, 0] = (projected_array[:, 0] -
                                 #(mouse_x/view_w))
        #projected_array[:, 1] = (projected_array[:, 1] -
                                 #(mouse_y/view_h))
        #distance_array = np.power(np.power(projected_array[:, 0], 2) +
                                  #np.power(projected_array[:, 1], 2), 0.5)
                                  
        
        m = self.projectionMatrix() * self.viewMatrix()
        
        print "\n\n\n", self.viewMatrix()
        print "\n", self.projectionMatrix()
        #quit()

        # Map 3D location to normalized device coordaintes: (-1,-1) and (1,1) are opposite corners of the view
        x=0; y=0; z=0
        pt = m.map(QtGui.QVector3D(self.mousePos.x(),self.mousePos.y(),0))
        proj_3d = pt*pts
        print len(pt*pts)    
        
        # Inverse mapping:
        pt = m.inverted()[0].map(QtGui.QVector3D(self.mousePos.x(),self.mousePos.y(),0))
        #proj_3d = pt*pts
        print len(pt*pts)                               
        
        for k in range(len(proj_3d)):
            for p in range(len(proj_3d[k])):
                plt.scatter(proj_3d[k][p].x(),proj_3d[k][p].y(), s=30)
            
        plt.show()
        
        
    def mousePos(self, ev):
        """ Store the position of the mouse press for later use.

        """
        self.mousePos = ev.pos()
        print "wtf.132.."

        
        
def plot_pyqt_3d(app, methods, files):
    global pts

    method = methods[0]

    w = MyGLViewWidget()
    #w = WorldView()
    w.setBackgroundColor([0,0,0])
    w.opts['distance'] = 20
    w.show()
    w.setWindowTitle('pyqtgraph example: GLScatterPlotItem')

    g = gl.GLGridItem()
    w.addItem(g)
        
    #Load data from file first
    ctr = 0
    start_rec = 5
    n_recs=1
    for file_ in files[start_rec:start_rec+n_recs]:
        file_temp = glob.glob(file_[:-4]+'_970*'+method+'.npy')
        coords=np.load(file_temp[0])[:20]
        
        scale_factor = 5
        coords=coords*scale_factor
        
        pts = []; size = []; color = []
        lines=[]; lines_colors=[]
        clr_ctr =0
        for i in range(len(coords)):
            locs = [coords[i][0]+scale_factor*(ctr%np.sqrt(n_recs)-np.sqrt(n_recs)/2.), 
                        coords[i][1]+scale_factor*(int(ctr/np.sqrt(n_recs))-np.sqrt(n_recs)/2.), 
                        coords[i][2]]
            pts.append(locs)
            size.append(1)
            color.append(pg.glColor((i,len(coords))))

            lines.append(locs)
            lines_colors.append(pg.glColor((i,len(coords))))
            
            clr_ctr+=1
            
        lines = np.vstack(lines)#.transpose()
        lines_colors = np.vstack(lines_colors)
        plt = gl.GLLinePlotItem(pos=lines, color=lines_colors, width=1., antialias=True)
        #w.addItem(plt)
            
        pts = np.array(pts)
        size = np.array(size)
        color = np.array(color)
        sp1 = gl.GLScatterPlotItem(pos=pts, size=size, color=color, pxMode=False)
        w.addItem(sp1)
        
        ctr+=1
    
    #print gl.itemsAt((1, 2, 3, 4))

    w.mousePos(w)
    print w.gluProject(w)
    x,y,z = gluProject(object.x,object.y,object.z,\
    glGetDoublev(GL_MODELVIEW_MATRIX),glGetDoublev(GL_PROJECTION_MATRIX),\
    glGetIntegerv(GL_VIEWPORT))

    ## Start Qt event loop unless running in interactive mode.
    QtGui.QApplication.instance().exec_()
    app.closeAllWindows()

def find_previous(array, value):
    temp = (np.abs(array-value)).argmin()
    if array[temp]>value: return temp-1
    else: return temp

def count_velocities(methods, files, main_dir):
    
    colors=['blue','red','green','black','cyan','pink','lightgreen','orange','yellow','brown','magenta']

    plot_counter=0
    for method in methods:
        counter=0
        for file_ in files:
            ax=plt.subplot(1,1,plot_counter+1)
            print method, plot_counter
            
            data_3d= np.load(main_dir + file_+'_'+method+'9701.npy')

            
            #ax = plt.subplot(2,2,plot_counter+1)
            #counter=0
            distances = []
            for k in range(len(data_3d)-1):
                #print data_3d[k]
                distances.append(scipy.spatial.distance.euclidean(data_3d[k], data_3d[k+1]))

            #bin_width = 
            x = np.arange(0,np.max(distances)*.3,.01)
            print x
            print len(distances)
            temp_hist = np.histogram(distances, bins = x)[0]
            
            plt.plot(x[:-1], temp_hist, color=colors[int(files[counter])-1], linewidth=2, label=files[counter])

            plt.legend(loc='upper right', fontsize=15)
            plt.title("Dim Red Method: "+method)
            
            counter+=1
        plot_counter+=1
    plt.suptitle("Velocity plots", fontsize=20)
    plt.show()
    
def count_pairwise_distances(methods, files, main_dir):

    colors=['blue','red','green','black','cyan','pink','lightgreen','orange','yellow','brown','magenta']

    method = methods[0]
    n_points = [2500, 5000, 7500, 10000, 5000, 10000, 10000]
    widths =   [2500, 2500, 2500,  2500, 5000,  5000, 10000]
    ctr=0
    for width, points in zip(widths, n_points):
        print width, points
        ax = plt.subplot(2,4,ctr+1)
        plot_counter=0
        for file_ in files[:10]:
            file_temp = glob.glob(file_[:-4]+'_970*'+method+'.npy')
            #if len(file_temp)==0: continue
            data_in=np.load(file_temp[0])
        
            print data_in.shape
            
            distances = scipy.spatial.distance.pdist(data_in[points-width:points],'euclidean')
            #print len(scipy.spatial.distance.pdist(data_in[0]))
        
            bin_width = .01
            x = np.arange(0,np.max(distances),bin_width)
            temp_hist = np.histogram(distances, bins = x)[0]

            plt.plot(x[:-1], temp_hist, color=colors[plot_counter%11], linewidth=2, label=file_temp[0][:-4].replace(main_dir,''))

            plot_counter+=1
            ax.set_xticks([])
            ax.set_yticks([])
            
        plt.legend(loc='upper right', fontsize=12)
        plt.title("Pairwise Dist Method\npts:"+str(points-width)+' to '+str(min(9700,points)))
        #break
        ctr+=1
        
    plt.show()


def count_pairwise_distances_all(methods, files, main_dir):

    colors=['blue','red','green','black','cyan','pink','lightgreen','orange','yellow','brown','magenta']

    method = methods[0]

    ax = plt.subplot(1,1,1)
    plot_counter=0
    ctr=0
    for file_ in files:
        file_temp = glob.glob(file_[:-4]+'_970*'+method+'.npy')
        if len(file_temp)==0: continue
        print "...loading: ", file_temp[0]
        distances = np.load(file_temp[0])[0:10000000]
   
        bin_width = .025
        max_dist = 10.0; #np.max(distances)
        x = np.arange(0,max_dist,bin_width)
        temp_hist = np.histogram(distances, bins = x)[0]

        plt.plot(x[:-1], temp_hist, color=colors[plot_counter%11], linewidth=2, label=file_temp[0][:-4].replace(main_dir,''), alpha=.8)
        plot_counter+=1

        #plt.scatter(0,np.max(temp_hist), s=100, color = 'black', alpha=.5)
        #plt.scatter(x[np.argmax(temp_hist)], 0, s=100, color = 'black', alpha=.5)

        ctr+=1

    plt.ylim(bottom=0)
    plt.xlim(left=0)        
    ax.set_xticks([])
    ax.set_yticks([])
    plt.legend(loc='upper right', fontsize=7)
    plt.title("Pairwise Dist Method - all Dimensions")
    
        
        
    plt.show()
   
    
def count_cube_distribution(methods, files, main_dir):
    
    plot_counter=0
    for method in methods:
        print method, plot_counter
        data_in=[]
        for file_ in files:
            data_in.append(np.load(main_dir + file_+'_'+method+'9701.npy'))
        
        colors=['blue','red','green','black','cyan','pink','lightgreen','orange','yellow','brown','magenta']
        
        grain = 100

        ax = plt.subplot(2,2,plot_counter+1)
        counter=0
        for data_3d in data_in:
            x = [np.max(data_3d[:,0]), np.min(data_3d[:,0])]
            y = [np.max(data_3d[:,1]), np.min(data_3d[:,1])]
            z = [np.max(data_3d[:,2]), np.min(data_3d[:,2])]
            
            grid=[]
            grid.append(np.linspace(x[0], x[1], grain))
            grid.append(np.linspace(y[0], y[1], grain))
            grid.append(np.linspace(z[0], z[1], grain))
            
            cube_3d = np.zeros((grain,grain,grain), dtype=np.int16)
            
            for k in range(len(data_3d)):
                coord = []
                for p in range(3):
                    coord.append(find_previous(grid[p], data_3d[k][p]))
                cube_3d[coord[0]][coord[1]][coord[2]]+=1

            print cube_3d.shape
            cube_3d_flat = cube_3d.reshape(grain*grain*grain)
            
            bin_width = 1
            x = np.arange(1,np.max(cube_3d_flat),bin_width)
            temp_hist = np.histogram(cube_3d_flat, bins = x)[0]
            
            plt.plot(x[:-1], temp_hist, color=colors[int(files[counter])-1], linewidth=2, label=files[counter])
            counter+=1

        plt.legend(loc='upper right', fontsize=15)
        plt.title("Dim Red Method: "+method+ " cube grain: "+str(grain))
        
        plot_counter+=1
    plt.xlim(1,10)
    #plt.xlim(0,10)
    plt.show()
    

def on_click(event):
    
    global coords, images_temp, ax, fig, cid
    
    n_pix = len(images_temp)
    
    if event.inaxes is not None:
        coords.append((event.ydata, event.xdata))
        for j in range(len(coords)):
            for k in range(3):
                for l in range(3):
                    images_temp[min(n_pix,int(coords[j][0])-1+k)][min(n_pix,int(coords[j][1])-1+l)]=np.max(images_temp)

        ax.imshow(images_temp)
        #plt.show()
        fig.canvas.draw()
                    #figManager = plt.get_current_fig_manager()
                    #figManager.window.showMaximized()
    else:
        print 'Exiting'
        plt.close()
        fig.canvas.mpl_disconnect(cid)

def Define_generic_mask_64(images_processed, main_dir):

    global coords, images_temp, ax, fig, cid
    
    images_temp = images_processed.copy()

    fig, ax = plt.subplots()

    if (os.path.exists(main_dir + 'genericmask_64.txt')==False):
        coords=[]

        ax.imshow(images_processed)#, vmin=0.0, vmax=0.02)
        ax.set_title("Compute generic (outside the brain) mask")
        #figManager = plt.get_current_fig_manager()
        #figManager.window.showMaximized()
        cid = fig.canvas.mpl_connect('button_press_event', on_click)
        plt.show()

        #******* MASK AND DISPLAY AREAS OUTSIDE GENERAL MASK 
        #Search points outside and black them out:
        all_points = []
        for i in range(len(images_processed[0])):
            for j in range(len(images_processed[0])):
                all_points.append([i,j])

        all_points = np.array(all_points)
        vertixes = np.array(coords) 
        vertixes_path = Path(vertixes)
        
        mask = vertixes_path.contains_points(all_points)
        counter=0
        coords_save=[]
        for i in range(len(images_processed[0])):
            for j in range(len(images_processed[0])):
                if mask[counter] == False:
                    images_processed[i][j]=0
                    coords_save.append([i,j])
                counter+=1

        fig, ax = plt.subplots()
        ax.imshow(images_processed)
        plt.show()
       
        genericmask_file = main_dir + 'genericmask_64.txt'
        np.savetxt(genericmask_file, coords_save)

        print "Finished Making General Mask"

    plt.close()

def mask_data(images_raw, main_dir, file_, data):
    print "Masking data..."
    img_rate = 150      #Frame rate of imaging
    n_pixels = len(data[0][0])      #number of pixels

    #Load General mask (removes background)
    #Define_generic_mask_64(data[0][500], main_dir)
    Define_generic_mask(images_raw[500], main_dir)

    generic_mask_file = []
    generic_mask_file = main_dir + 'genericmask.txt'
    #generic_mask_file = main_dir + 'genericmask_64.txt'
    if (os.path.exists(generic_mask_file)==True):
        generic_coords = np.loadtxt(generic_mask_file)
    else:
        #generic_coords = Define_generic_mask_64(images_filtered, main_dir, file_, img_rate, n_pixels)
        generic_coords = Define_generic_mask(images_filtered, main_dir, file_, img_rate, n_pixels)

    generic_mask_indexes=np.zeros((n_pixels,n_pixels))
    for i in range(len(generic_coords)):
        generic_mask_indexes[int(generic_coords[i][0])][int(generic_coords[i][1])] = True

    for k in range(len(data)):
        temp_array = []
        for i in range(0, len(data[k]),1):
            temp = np.ma.array(data[k][i], mask=generic_mask_indexes, fill_value = 0., hard_mask = True)
            temp_array.append(temp)
        data[k] = np.ma.array(temp_array)

def Define_generic_mask(images_processed, main_dir):

    global coords, images_temp, ax, fig, cid
    
    images_temp = images_processed.copy()

    fig, ax = plt.subplots()

    if (os.path.exists(main_dir + 'genericmask.txt')==False):
        coords=[]

        ax.imshow(images_processed)#, vmin=0.0, vmax=0.02)
        ax.set_title("Compute generic (outside the brain) mask")
        #figManager = plt.get_current_fig_manager()
        #figManager.window.showMaximized()
        cid = fig.canvas.mpl_connect('button_press_event', on_click)
        plt.show()

        #******* MASK AND DISPLAY AREAS OUTSIDE GENERAL MASK 
        #Search points outside and black them out:
        all_points = []
        for i in range(len(images_processed[0])):
            for j in range(len(images_processed[0])):
                all_points.append([i,j])

        all_points = np.array(all_points)
        vertixes = np.array(coords) 
        vertixes_path = Path(vertixes)
        
        mask = vertixes_path.contains_points(all_points)
        counter=0
        coords_save=[]
        for i in range(len(images_processed[0])):
            for j in range(len(images_processed[0])):
                if mask[counter] == False:
                    images_processed[i][j]=0
                    coords_save.append([i,j])
                counter+=1

        fig, ax = plt.subplots()
        ax.imshow(images_processed)
        plt.show()
       
        genericmask_file = main_dir + 'genericmask.txt'
        np.savetxt(genericmask_file, coords_save)

        print "Finished Making General Mask"

    plt.close()


def make_vids(data, file_):
    
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=50, metadata=dict(artist='Me'), bitrate=6400)

    fig = plt.figure() # make figure

    im = []
    for k in range(len(data)):
        ax = plt.subplot(1,2,k+1)

        ax.get_xaxis().set_visible(False)
        ax.yaxis.set_ticks([])
        ax.yaxis.labelpad = 0
        Ev_max = np.max(np.abs(data[k]))*.9
        #v_min = -v_max 
        
        im.append([])
        im[k] = plt.imshow(data[k][0], cmap=plt.get_cmap('gray'), vmin = -Ev_max, vmax = Ev_max, interpolation='none')#, vmin=0, vmax=v_max)

    #function to update figure
    def updatefig(j):
        print "... frame: ", j
        plt.suptitle("Frame: "+str(j)+"\n"+str(round(float(j)/150,2))+"sec")

        # set the data in the axesimage object
        for k in range(len(data)):
            im[k].set_array(data[k][j])

        # return the artists set
        return im
        
    ani = animation.FuncAnimation(fig, updatefig, frames=range(len(data[0])), interval=10, blit=False, repeat=True)
    ani.save(file_+'.mp4', writer=writer)
    plt.show()

def process_data(files, main_dir, dim_red_methods):

    for file_ in files: #[len(files)/2:]:
        file_ = file_[:-4]
        print "Processing file: ", file_
        
        #Read .tif
        if (os.path.exists(file_+'.npy')==False):
            images_raw = tiff.imread(file_+'.tif')
            np.save(file_, images_raw)

        #Filter .tif
        #if True:
        if (os.path.exists(file_+'_filtered.npy')==False):
            images_raw = np.load(file_+'.npy')

            print "Temporal filtering ..."
            from scipy.signal import butter, lfilter, filtfilt

            def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
                b, a = butter_bandpass(lowcut, highcut, fs, order=order)
                y = filtfilt(b, a, data,axis=0)
                return y

            def butter_bandpass(lowcut, highcut, fs, order=5):
                nyq = 0.5 * fs
                low = lowcut / nyq
                high = highcut / nyq
                b, a = butter(order, [low, high], btype='band')
                return b, a

            lowcut = .5
            highcut=6
            img_rate = 150        #Frame rate of imaging

            print "...band pass filtering: ", lowcut, "hz to ", highcut, "hz"
            images_filtered = np.array(butter_bandpass_filter(np.float32(images_raw), lowcut, highcut, img_rate, order = 2))

            np.save(file_+'_filtered', np.float32(images_filtered))

        images_raw = np.load(file_+'.npy')
        images_filtered = np.load(file_+'_filtered.npy')
        print images_filtered.shape

        print "...shifting DFF ..."
        print np.nanmin(images_filtered)

        #images_filtered3 = images_filtered + 10*images_raw[300]
        #images_filtered4 = images_filtered - np.ndarray.min(images_filtered, axis=0)


        #DFF Computation
        data = []
        print "...computing DFF..."
        print "...cat's method..."
        images_filtered1 = images_filtered - np.nanmin(images_filtered) #+ images_raw[300]
        baseline = np.average(images_filtered1[300:], axis=0)       #Cat's method: uses minimum offset for processing
        data1 = (images_filtered1[300:]-baseline)/baseline
        data.append(data1[:])
        #data.append(images_filtered1[300:])

        #print "...Allen's method..."
        #images_filtered2 = images_filtered + images_raw[300]
        #baseline = np.nanmean(images_filtered2[300:], axis=0)       #Allen/Matthieu method: uses raw frame[300] to offset data.
        #data2 = (images_filtered2[300:]-baseline)/baseline
        #data.append(data2[:])

        #print "...Jeff's's method..."
        #images_filtered3 = images_filtered + np.average(images_raw) #NEED TO DOUBLE CHECK THIS
        #baseline = np.average(images_filtered3[300:], axis=0)
        #data3 = (images_filtered3[300:]-baseline)/baseline
        #data.append(data3[:])

        ##baseline = np.average(images_filtered4[300:], axis=0)
        ##data4 = (images_filtered4[300:]-baseline)/baseline
        ##data.append(data4[:])
        #data.append(images_filtered4[300:])


        #******* MASK DATA ********
        mask_data(images_raw, main_dir, file_, data)
            
        #***********GENERATE ANIMATIONS
        if False: make_vids(make_vids(data, file_))


        #data[0] = data[0][:100]
        #********* DIMENSIONALITY REDUCTION **************
        #smooth data to 64 x 64 pixels
        from skimage.measure import block_reduce
        print data[0].shape
        data[0] = block_reduce(data[0], block_size=(1,2,2), func=np.mean)
        print data[0].shape

        #Vectorize data
        temp_data = data[0].reshape(data[0].shape[0],4096)

        #dim_red_methods = [0,1,2,3]
        dim_red_methods = [1]
        for dim_red_method in dim_red_methods:
            dim_reduction(temp_data, dim_red_method, file_)      #Methods: 0: MDS; 1:tSNE; 2: PCA; 3: Barnes-Hut tSNE

