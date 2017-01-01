"""PyQt OpenGL example, modified from PySide example, which was a port of the opengl/hellogl
example from Qt v4.x"""

from __future__ import division

#from cell_gui import *

from scipy.stats import binned_statistic
from operator import truediv, sub
import struct, array, csv
import h5py
import glob
from scipy.interpolate import interp1d
from scipy.signal import butter, lfilter

import json, math
from pprint import pprint
from scipy.io import mmread

import sys, time, os
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import show, plot
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib.colors import LinearSegmentedColormap
from pylab import rcParams

from OpenGL.GL import *
from OpenGL.GLU import *
#from OpenGL.WGL import *

#import IPython

# instantiate an IPython embedded shell which shows up in the terminal on demand
# and on every exception:
#from IPython.terminal.ipapp import load_default_config
#from IPython.terminal.embed import InteractiveShellEmbed
#config = load_default_config()
# automatically call the pdb debugger after every exception, override default config:
#config.TerminalInteractiveShell.pdb = True
#ipshell = InteractiveShellEmbed(display_banner=False, config=config)

from PyQt4 import QtCore, QtGui, QtOpenGL, uic
from PyQt4.QtCore import Qt

MainUi, MainUiBase = uic.loadUiType('main.ui')

from OpenGL import GL, GLU
import numpy as np

'''
def normdeg(angle):
    return angle % 360
'''

brb   = {'red':   ((0.0, 0.0, 0.0),
                   (0.5, 0.1, 0.1),
                   (1.0, 1.0, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 1.0, 1.0),
                   (0.5, 0.1, 0.1),
                   (1.0, 0.0, 0.0))
        }

custom_jet =   {'red':   ((0., 0, 0), (0.35, 0, 0), (0.66, 1, 1), (0.89,1, 1),
                         (1, 0.5, 0.5)),
               'green': ((0., 0, 0), (0.125,0, 0), (0.375,0.5, .5), (0.64,1.0, 1.0),
                         (0.91,1.0,1.0), (1, 0, 0)),
               'blue':  ((0., 0.5, 0.5), (0.11, 1, 1), (0.34, 1, 1), (0.65,0, 0),
                         (1, 0, 0))}

RED = 255, 0, 0
ORANGE = 255, 127, 0
YELLOW = 255, 255, 0
GREEN = 0, 255, 0
CYAN = 0, 255, 255
LIGHTBLUE = 0, 127, 255
BLUE = 0, 0, 255
VIOLET = 127, 0, 255
MAGENTA = 255, 0, 255  
GREY = 85, 85, 85
WHITE = 255, 255, 255
DARK_GREY = 30, 30, 30
BLACK = 0,0,0
CMAP = np.array([RED, ORANGE, YELLOW, GREEN, CYAN, LIGHTBLUE, BLUE, VIOLET, MAGENTA,
                 GREY, WHITE, DARK_GREY, BLACK], dtype=np.uint8)

COLOUR_DICT = {'pyramid1': 1, 'pyramid2': 2, 'pyramid3': 3, 'pyramid4': 4, 'pyramid21': 5,
               'pyramid26': 6, 'pyramid28': 7, 'bc1': 0, 'bc2': 8, 'bc3': 1, 'bc4': 9, 'bc5': 10, 'bc6': 11,
	       'pyramid23': 1, 'pyramid5': 2, 'pyramid6': 3, 
	       'basket23': 5, 'basket4': 10, 'basket5': 7, 'basket6': 8}

RATCELL_DICT = {'pyramid1': 0, 'pyramid2': 1, 'pyramid3': 2, 'pyramid4': 3, 'pyramid21': 4,
               'pyramid26': 5, 'pyramid28': 6, 'bc1': 7, 'bc2': 8, 'bc3': 0, 'bc4': 9, 'bc5': 10, 'bc6': 11}

RATPOP_NAMES = ['pyramid1', 'pyramid2', 'pyramid3', 'pyramid4', 'pyramid21',
               'pyramid26', 'pyramid28', 'bc1', 'bc2', 'bc4', 'bc5', 'bc6']



MOUSECELL_DICT = {'pyramid23': 0, 'pyramid4': 1, 'pyramid5': 2, 'pyramid6': 3,
                 'basket23': 4, 'basket4': 5, 'basket5': 6, 'basket6': 7}

MOUSEPOP_NAMES = ['pyramid23', 'pyramid4', 'pyramid5', 'pyramid6',
                 'basket23', 'basket4', 'basket5', 'basket6']



class MainWindow(QtGui.QMainWindow):
    
    def __init__(self, parent=None):
        
        QtGui.QWidget.__init__(self, parent)
        #Button menu
        self.ui = MainUi()
        self.ui.setupUi(self) # lay it out
        self.glwindow = GLWindow(parent=None)

    @QtCore.pyqtSlot()
    def on_plotSingleTrialEphys_clicked(self):
        self.glwindow.show()
        
        print "plotFixedPts"

        glwin = self.glwindow
        
        glwin.glWidget.plot_soma = 1
        glwin.glWidget.plot_vertices = 1
        glwin.glWidget.plot_frame=1
        #glwin.glWidget.file_name = self.file_name

        primitives = draw_primitives()
        size = 2
        cmap = mpl.cm.get_cmap('jet')

        colors_temp = []
        vertices_locs = []

        #Use single dim red call - separate data out manually
        print len(self.single_trial)
        primitives.spheres(self.single_trial, size) #Generate spheres at 3D point locations

        glwin.glWidget.quad_points = np.vstack([primitives.points_quads]) 
        glwin.glWidget.quad_colors = np.vstack([primitives.points_colors])

        print "... done..."
        #For single movie repeats - connect Nodes sequentially
        print "Generating vertices..."
        for k in range(len(self.single_trial)-1):
            vertices_locs.append(self.single_trial[k])
            vertices_locs.append(self.single_trial[k+1])

            vertex_color = cmap(float(k)/len(self.single_trial))
            
            colors_temp.append([vertex_color[0],vertex_color[1],vertex_color[2],1.0])
            colors_temp.append([vertex_color[0],vertex_color[1],vertex_color[2],1.0])

        #print vertices_locs

        glwin.glWidget.vertices_points = np.vstack([vertices_locs])
        glwin.glWidget.vertices_colors = np.vstack([colors_temp])
    
        print "...done..."
        
        #****** PLOT FRAME ******
        frame_scale=.75
        primitives.simFrame(frame_scale)
        
        glwin.glWidget.points_frame = np.vstack([primitives.points_frame])
        glwin.glWidget.colours_frame  = np.vstack([primitives.colors_frame])
        
                
        glwin.glWidget.updateGL()
        
    @QtCore.pyqtSlot()
    def on_plotMultiTrialEphys_clicked(self):
        self.glwindow.show()
        
        print "MultiTrialEphys"

        glwin = self.glwindow
        
        glwin.glWidget.plot_soma = 1
        glwin.glWidget.plot_vertices = 1
        glwin.glWidget.plot_frame=1
        #glwin.glWidget.file_name = self.file_name

        primitives = draw_primitives()
        size = 2
        cmap = mpl.cm.get_cmap('jet')

        colors_temp = []
        vertices_locs = []

        ##Plot movie repeat same time vector distances/locations
        print "multi_trial shape in (pre-chunk): ", np.array(self.multi_trial).shape

        #Plot Nodes
        for r in range(self.repeats):
            primitives.spheres(self.multi_trial[r],size) #Generate spheres at 3D point locations

            if r>0:
                glwin.glWidget.quad_points = np.vstack((glwin.glWidget.quad_points, primitives.points_quads)) 
                glwin.glWidget.quad_colors = np.vstack((glwin.glWidget.quad_colors, primitives.points_colors ))
            else:
                glwin.glWidget.quad_points = np.vstack([primitives.points_quads]) 
                glwin.glWidget.quad_colors = np.vstack([primitives.points_colors])

        #Plot Vertices
        for r in range(self.repeats-1):
            print "Vertices trial: ", r
            vertices_locs=[]
            for k in range(len(self.multi_trial[r])):
                vertices_locs.append(self.multi_trial[r][k])
                vertices_locs.append(self.multi_trial[r+1][k])

                vertex_color = cmap(float(k)/len(self.multi_trial[r]))
                
                colors_temp.append([vertex_color[0],vertex_color[1],vertex_color[2],1.0])
                colors_temp.append([vertex_color[0],vertex_color[1],vertex_color[2],1.0])
            
            #print vertices_locs

            if r>0:
                print glwin.glWidget.vertices_points.shape, np.array(vertices_locs).shape
                glwin.glWidget.vertices_points = np.vstack((glwin.glWidget.vertices_points, vertices_locs ))
                glwin.glWidget.vertices_colors = np.vstack((glwin.glWidget.vertices_colors, colors_temp ))
            else:
                glwin.glWidget.vertices_points = np.vstack([vertices_locs])
                glwin.glWidget.vertices_colors = np.vstack([colors_temp])

                #colors_temp.append([vertex_color[0],vertex_color[1],vertex_color[2],1.0])
                #colors_temp.append([vertex_color[0],vertex_color[1],vertex_color[2],1.0])
            
            #glwin.glWidget.vertices_points = np.vstack([vertices_locs])
            #glwin.glWidget.vertices_colors = np.vstack([colors_temp])
        
        
        print "...done..."
        #****** PLOT FRAME ******
        frame_scale=.75
        primitives.simFrame(frame_scale)
        
        glwin.glWidget.points_frame = np.vstack([primitives.points_frame])
        glwin.glWidget.colours_frame  = np.vstack([primitives.colors_frame])
        
        
        glwin.glWidget.updateGL()
        

    @QtCore.pyqtSlot()
    def on_plotSingleTrialImaging_clicked(self):
        self.glwindow.show()
        
        print "plotFixedPts"

        glwin = self.glwindow
        
        glwin.glWidget.plot_soma = 1
        glwin.glWidget.plot_vertices = 1
        glwin.glWidget.plot_frame=1

        primitives = draw_primitives()
        size = 2
        cmap = mpl.cm.get_cmap('jet')

        colors_temp = []
        vertices_locs = []
        offset = 100.

        #Plot sphere/cube locations
        print len(self.single_trial), self.time_segs
        primitives.spheres(self.single_trial,size) #Generate spheres at 3D point locations

        glwin.glWidget.quad_points = np.vstack([np.array(primitives.points_quads)]) 
        glwin.glWidget.quad_colors = np.vstack([primitives.points_colors])
           
        #Plot vertices
        print "Generating vertices..."
        for k in range(len(self.single_trial)-1):
            vertices_locs.append(self.single_trial[k])
            vertices_locs.append(self.single_trial[k+1])

            vertex_color = cmap(float(k)/len(self.single_trial))
            
            colors_temp.append([vertex_color[0],vertex_color[1],vertex_color[2],1.0])
            colors_temp.append([vertex_color[0],vertex_color[1],vertex_color[2],1.0])
        
        glwin.glWidget.vertices_points = np.vstack([np.array(vertices_locs)])
        glwin.glWidget.vertices_colors = np.vstack([colors_temp])
    
        print "...done..."

        #****** PLOT FRAME ******
        frame_scale=.75
        primitives.simFrame(frame_scale)
        
        glwin.glWidget.points_frame = np.vstack([primitives.points_frame])
        glwin.glWidget.colours_frame  = np.vstack([primitives.colors_frame])
                
        glwin.glWidget.updateGL()
        
    @QtCore.pyqtSlot()
    def on_plotMultiTrialImaging_clicked(self):
        self.glwindow.show()
        
        print "plotFixedPts"

        glwin = self.glwindow
        
        glwin.glWidget.plot_soma = 1
        glwin.glWidget.plot_vertices = 1
        glwin.glWidget.plot_frame=1

        primitives = draw_primitives()
        size = 2
        cmap = mpl.cm.get_cmap('jet')

        colors_temp = []
        vertices_locs = []
        offset = 100.
        
        #Plot sphere/cube locations
        
        from rigid_transform import rigid_transform_3D

        trials = []
        counter = 0
        for file_name in (self.file_names):

            print "Loading: ", file_name
            trial = np.load(self.work_dir + file_name+'_nodes_vertices_method_'+str(self.method)+'.npy')
            #if counter>0:
            #    trial = trial -trials[0][0]

            #Align to 1st point
            if counter>0:
                R, t = rigid_transform_3D(trial.copy(),trials[0])
                print "Rotation: ", R
                print "Transformation: ", t
                trial = np.dot(R, trial.T).T 

            trials.append(trial)
                
            counter+=1
        #quit()
        
        counter=0
        for trial in trials:
            primitives.spheres(trial,size) #Generate spheres at 3D point locations
            if counter>0:
                glwin.glWidget.quad_points = np.vstack((glwin.glWidget.quad_points, np.array(primitives.points_quads))) 
                glwin.glWidget.quad_colors = np.vstack((glwin.glWidget.quad_colors, primitives.points_colors ))
            else:
                glwin.glWidget.quad_points = np.vstack([np.array(primitives.points_quads)]) 
                glwin.glWidget.quad_colors = np.vstack([primitives.points_colors])
           
            print "Generating vertices..."
            for k in range(len(trial)-1):
                vertices_locs.append(trial[k])
                vertices_locs.append(trial[k+1])


                vertex_color = cmap(float(k)/len(trial))
                
                colors_temp.append([vertex_color[0],vertex_color[1],vertex_color[2],1.0])
                colors_temp.append([vertex_color[0],vertex_color[1],vertex_color[2],1.0])
            
            if counter>0:
                glwin.glWidget.vertices_points = np.vstack((glwin.glWidget.vertices_points, np.array(vertices_locs)))
                glwin.glWidget.vertices_colors = np.vstack((glwin.glWidget.vertices_colors, colors_temp ))
            else:
                glwin.glWidget.vertices_points = np.vstack([np.array(vertices_locs)])
                glwin.glWidget.vertices_colors = np.vstack([colors_temp])
        
            print "...done..."
        
            counter+=1

        #****** PLOT FRAME ******
        frame_scale=.75
        primitives.simFrame(frame_scale)
        
        glwin.glWidget.points_frame = np.vstack([primitives.points_frame])
        glwin.glWidget.colours_frame  = np.vstack([primitives.colors_frame])
                
        glwin.glWidget.updateGL()
        

    @QtCore.pyqtSlot()
    def on_mousePick_clicked(self):
        print "mousePick"

        self.glwindow.show()

        glwin = self.glwindow
        
        glwin.glWidget.plot_soma = 1
        glwin.glWidget.pick()
                
        #glwin.glWidget.updateGL()
        
#***************************************************************************************
#************************************ GLWindow Class ***********************************
#***************************************************************************************
                       
class GLWindow(QtGui.QWidget):
    
    def __init__(self, parent=None):
        
        QtGui.QWidget.__init__(self, parent)
        self.glWidget = GLWidget(parent=self)
        mainLayout = QtGui.QHBoxLayout()
        mainLayout.addWidget(self.glWidget)
        self.setLayout(mainLayout)
        self.setWindowTitle(self.tr("OpenGL test"))


#***************************************************************************************
#************************************ GLWidget Class ***********************************
#***************************************************************************************

class GLWidget(QtOpenGL.QGLWidget):
    
    def __init__(self, parent=None):
        
        QtOpenGL.QGLWidget.__init__(self, parent)
        self.lastPos = QtCore.QPoint()
        self.focus = np.float32([0, 0, 0]) # init camera focus
        self.axes = 'both' # display both mini and focal xyz axes by default

        format = QtOpenGL.QGLFormat()
        #format.setVersion(3, 0) # not available in PyQt 4.7.4
        # set to color index mode, unsupported in OpenGL >= 3.1, don't know how to load
        # GL_ARB_compatibility extension, and for now, can't force OpenGL 3.0 mode.
        # Gives "QGLContext::makeCurrent(): Cannot make invalid context current." error:
        #format.setRgba(False)
        
        format.setDoubleBuffer(True) # works fine
        self.setFormat(format)
        #QtOpenGL.QGLFormat.setDefaultFormat(format)
        
        '''
        c = QtGui.qRgb
        cmap = [c(255, 0, 0), c(0, 255, 0), c(0, 0, 255), c(255, 255, 0), c(255, 0, 255)]
        colormap = QtOpenGL.QGLColormap()
        colormap.setEntries(cmap)
        self.setColormap(colormap)
        '''
        self.npoints = 100000
        if self.npoints > 2**24-2: # the last one is the full white bg used as a no hit
            raise OverflowError("Can't pick from more than 2**24-2 sids")
        self.points = np.float32(np.random.random((self.npoints, 3))) - 0.5
        self.sids = np.arange(self.npoints)
        self.nids = self.sids % len(CMAP)
        self.colors = CMAP[self.nids] # uint8
        # encode sids in RGB
        r = self.sids // 256**2
        rem = self.sids % 256**2 # remainder
        g = rem // 256
        b = rem % 256
        self.rgbsids = np.zeros((self.npoints, 3), dtype=np.uint8)
        self.rgbsids[:, 0] = r
        self.rgbsids[:, 1] = g
        self.rgbsids[:, 2] = b
        #print self.rgbsids
        #print self.rgbsids.dtype



    def minimumSizeHint(self):
        return QtCore.QSize(1000, 800)

    def sizeHint(self):
        return QtCore.QSize(1000, 800)

    def initializeGL(self):
        bkgr = 1
        if bkgr==0:
            GL.glClearColor(0.0, 0.0, 0.0, 1.0) # same as default
        else:
            GL.glClearColor(1.0, 1.0, 1.0, 1.0)

        GL.glClearDepth(10.0) # same as default
        GL.glEnable(GL.GL_DEPTH_TEST) # display points according to occlusion, not order of plotting
        #GL.glEnable(GL.GL_POINT_SMOOTH) # doesn't seem to work right, proper way to antialiase?
        #GL.glEnable(GL.GL_LINE_SMOOTH) # works better
        #GL.glPointSize(1.5) # truncs to the nearest pixel if antialiasing is off
        GL.glShadeModel(GL.GL_FLAT)
        #GL.glEnable(GL.GL_CULL_FACE) # only useful for solids
        GL.glTranslate(0, 750, -3000) # init camera distance from origin

    def paintGL(self):
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)

        GL.glEnable(GL_BLEND)
        GL.glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

        # Don't load identity matrix. Do all transforms in place against current matrix
        # and take advantage of OpenGL's state-machineness.
        # Sure, you might get round-off error over time, but who cares? If you wanna return
        # to a specific focal point on 'f', that's when you really need to first load the
        # identity matrix before doing the transforms
        #GL.glLoadIdentity() # loads identity matrix into top of matrix stack

        # viewing transform for camera: where placed, where pointed, which way is up:
        #GLU.gluLookAt()
        #GL.glScale() # modelling transformation, lets you stretch your objects

        GL.glEnableClientState(GL.GL_COLOR_ARRAY);
        GL.glEnableClientState(GL.GL_VERTEX_ARRAY);
           

        #PLOT NODES
        if self.plot_soma==1:
            ##GL.glColorPointerub(self.triangle_colours) # unsigned byte, ie uint8
            ##GL.glVertexPointerf(self.triangle_points) # float32
            ##GL.glDrawArrays(GL.GL_TRIANGLES, 0, len(self.triangle_points)*3)

            GL.glColorPointer(4, GL.GL_FLOAT, 0, self.quad_colors) # float color 0 .. 1.0
            GL.glVertexPointerf(self.quad_points) # float32
            GL.glDrawArrays(GL.GL_QUADS, 0, len(self.quad_points)*4)

        #PLOT VERTICES
        if self.plot_vertices ==1:
            #GL.glColorPointerub(self.vertices_colors) # unsigned byte, ie uint8
            GL.glColorPointer(4, GL.GL_FLOAT, 0, self.vertices_colors) # float color 0 .. 1.0
            GL.glVertexPointerf(self.vertices_points) # float32
            GL.glDrawArrays(GL.GL_LINES, 0, len(self.vertices_points))

        #PLOT TITLE/TEXT
        if False:
            GL.glColor3ub(0, 0, 0)
            self.renderText (len(self.file_name),0,1000, self.file_name)

        #PLOT FRAME
        if self.plot_frame==1:
            GL.glColorPointerub(self.colours_frame) # unsigned byte, ie uint8
            GL.glVertexPointerf(self.points_frame) # float32
            GL.glDrawArrays(GL.GL_LINES, 0, len(self.points_frame))

        if self.axes: # paint xyz axes
            GL.glClear(GL.GL_DEPTH_BUFFER_BIT) # make axes paint on top of data points
            #if self.axes in ['both', 'mini']:
            self.paint_mini_axes()
            #if self.axes in ['both', 'focal']:
            self.paint_focal_axes()

        # might consider using buffer objects for even more speed (less unnecessary vertex
        # data from ram to vram, I think). Apparently, buffer objects don't work with
        # color arrays?

        #GL.glFlush() # forces drawing to begin, only makes difference for client-server?
        self.swapBuffers() # doesn't seem to be necessary, even though I'm in double-buffered
                           # mode with the back buffer for RGB sid encoding, but do it anyway
                           # for completeness

    def resizeGL(self, width, height):
        GL.glViewport(0, 0, width, height)
        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glLoadIdentity()
        # fov (deg) controls amount of perspective, and as a side effect initial apparent size
        GLU.gluPerspective(45, width/height, 0.01, 1000000) # fov, aspect, nearz & farz
                                                           # clip planes
        GL.glMatrixMode(GL.GL_MODELVIEW)
    

    def paint_mini_axes(self):
        """Paint mini xyz axes in bottom left of widget"""
        w, h = self.width(), self.height()
        vt = self.getTranslation() # this is in eye coordinates
        GL.glViewport(0, 0, w//4, h//4) # mini viewport at bottom left of widget
        self.setTranslation((-1, -.5, -3)) # draw in center of this mini viewport
        self.paint_axes()
        self.setTranslation(vt) # restore translation vector to MV matrix
        GL.glViewport(0, 0, w, h) # restore full viewport

    def paint_focal_axes(self):
        """Paint xyz axes proportional in size to sigma, at focus"""
        GL.glTranslate(*self.focus) # translate to focus
        #self.paint_axes(self.sigma)
        GL.glTranslate(*-self.focus) # translate back

    def update_focal_axes(self):
        """Called every time sigma is changed in main spyke window"""
        #self.update_sigma()
        self.updateGL()

    def paint_axes(self, l=1):
        """Paint axes at origin, with lines of length l"""
        GL.glBegin(GL.GL_LINES)
        GL.glColor3f(1, 0, 0) # red x axis
        GL.glVertex3f(0, 0, 0)
        GL.glVertex3f(l, 0, 0)
        GL.glColor3f(0, 1, 0) # green y axis
        GL.glVertex3f(0, 0, 0)
        GL.glVertex3f(0, l, 0)
        GL.glColor3f(0, 0, 1) # blue z axis
        GL.glVertex3f(0, 0, 0)
        GL.glVertex3f(0, 0, l)
        GL.glEnd()

    def get_MV(self):
        """Return modelview matrix"""
        return GL.glGetDoublev(GL.GL_MODELVIEW_MATRIX) # I think this acts like a copy

    def set_MV(self, MV):
        GL.glLoadMatrixd(MV)

    MV = property(get_MV, set_MV)

    # modelview matrix is column major, so we work on columns instead of rows
    def getViewRight(self):
        """View right vector: 1st col of modelview matrix"""
        return self.MV[:3, 0]

    def getViewUp(self):
        """View up vector: 2nd col of modelview matrix"""
        return self.MV[:3, 1]

    def getViewNormal(self):
        """View normal vector: 3rd col of modelview matrix"""
        return self.MV[:3, 2]

    def getTranslation(self):
        """Translation vector: 4th row of modelview matrix"""
        return self.MV[3, :3]

    def setTranslation(self, vt):
        """Translation vector: 4th row of modelview matrix"""
        MV = self.MV
        MV[3, :3] = vt
        self.MV = MV

    def getDistance(self):
        v = self.getTranslation()
        #return np.sqrt((v**2).sum()) # from data origin
        return np.sqrt(((v-self.focus)**2).sum()) # from focus

    def pan(self, dx, dy):
        """Translate along view right and view up vectors"""
        d = self.getDistance()
        vr = self.getViewRight()
        vr *= dx*d
        GL.glTranslate(vr[0], vr[1], vr[2])
        vu = self.getViewUp()
        vu *= dy*d
        GL.glTranslate(vu[0], vu[1], vu[2])

    def zoom(self, dr):
        """Translate along view normal vector"""
        d = self.getDistance()
        vn = self.getViewNormal()
        vn *= dr*d
        GL.glTranslate(vn[0], vn[1], vn[2])

    def pitch(self, dangle): # aka elevation
        """Rotate around view right vector"""
        vr = self.getViewRight()
        GL.glTranslate(*self.focus)
        GL.glRotate(dangle, *vr)
        GL.glTranslate(*-self.focus)

    def yaw(self, dangle): # aka azimuth
        """Rotate around view up vector"""
        vu = self.getViewUp()
        GL.glTranslate(*self.focus)
        GL.glRotate(dangle, *vu)
        GL.glTranslate(*-self.focus)

    def roll(self, dangle):
        """Rotate around view normal vector"""
        vn = self.getViewNormal()
        GL.glTranslate(*self.focus)
        GL.glRotate(dangle, *vn)
        GL.glTranslate(*-self.focus)

    def panTo(self, p=None):
        """Translate along view right and view up vectors such that data point p is
        centered in the viewport. Not entirely sure why or how this works, figured
        it out using guess and test"""
        if p == None:
            p = self.focus
        MV = self.MV
        vr = self.getViewRight()
        vu = self.getViewUp()
        p = -p
        x = np.dot(p, vr) # dot product
        y = np.dot(p, vu)
        MV[3, :2] = x, y # set first two entries of 4th row to x, y
        self.MV = MV

    def pick(self):
        globalPos = QtGui.QCursor.pos()
        pos = self.mapFromGlobal(globalPos)
        width = self.size().width()
        height = self.size().height()
        x = pos.x()
        y = height - pos.y()
        print x, y
        if not (0 <= x < width and 0 <= y < height):
            print('cursor out of range')
            return
        '''
        # for speed, 1st check if there are any non-black pixels around cursor:
        GL.glReadBuffer(GL.GL_FRONT)
        frontbuffer = GL.glReadPixelsub(0, 0, width, height, GL.GL_RGB) # unsigned byte
        #rgb = frontbuffer[y-1:y+2, x-1:x+2] # +/- 1 pix
        rgb = frontbuffer[y, x]
        #print('frontbuffer:')
        #print rgb
        if (rgb == 0).all():
            print('nothing to return')
            return # nothing to pick
        '''
        
        # drawing encoded RGB values to back buffer
        #GL.glDrawBuffer(GL_BACK) # shouldn't be necessary, defaults to back
        GL.glClearColor(1.0, 1.0, 1.0, 1.0) # highest possible RGB means no hit
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)

        GL.glEnableClientState(GL.GL_COLOR_ARRAY);
        GL.glEnableClientState(GL.GL_VERTEX_ARRAY);
        GL.glColorPointerub(self.rgbsids) # unsigned byte, ie uint8
        GL.glVertexPointerf(self.points) # float32

        GL.glDrawArrays(GL.GL_POINTS, 0, self.npoints) # to back buffer

        GL.glClearColor(0.0, 0.0, 0.0, 1.0) # restore to default black

        # grab back buffer
        GL.glReadBuffer(GL.GL_BACK)

        # find rgb at cursor coords, decode sid
        # unsigned byte, x, y is bottom left:
        print x, y
        print GL.GL_RGB
        backbuffer = GL.glReadPixelsub(x, y, 1, 1, GL.GL_RGB)

        #rgb = backbuffer[0, 0]
        #r, g, b = rgb
        r,g,b = 255,255,255
        
        sid = r*256**2 + g*256 + b
        if sid == 2**24 - 1:
            #print('no hit')
            return
        nid = sid % len(CMAP)
        color = CMAP[nid]
        print('backbuffer.shape: %r' % (backbuffer.shape,))
        print rgb
        print('sid, nid, color: %d, %d, %r' % (sid, nid, color))
        '''
        # TODO: this isn't even necessary, since back buffer gets overdrawn anyway on next
        # updateGL():
        # restore front buffer to back
        GL.glReadBuffer(GL.GL_FRONT)
        #GL.glDrawBuffer(GL_BACK) # shouldn't be necessary, defaults to back
        GL.glRasterPos(0, 0) # destination position in draw buffer
        GL.glCopyPixels(0, 0, width, height, GL.GL_COLOR) # copies from read to draw buffer
        '''
        #self.swapBuffers() # don't even need to swap buffers, cuz we haven't changed the scene
        return sid

 
    def mousePressEvent(self, event):
        self.lastPos = QtCore.QPoint(event.pos())

    def mouseMoveEvent(self, event):
        buttons = event.buttons()
        modifiers = event.modifiers()
        dx = event.x() - self.lastPos.x()
        dy = event.y() - self.lastPos.y()

        if buttons == QtCore.Qt.LeftButton:
            if modifiers == Qt.ControlModifier:
                self.roll(-0.5*dx - 0.5*dy)
            elif modifiers == Qt.ShiftModifier:
                self.pan(dx/600, -dy/600) # qt viewport y axis points down
            else:
                self.yaw(0.5*dx)
                self.pitch(0.5*dy)
        elif buttons == QtCore.Qt.RightButton:
            self.zoom(-dy/500) # qt viewport y axis points down

        self.updateGL()
        self.lastPos = QtCore.QPoint(event.pos())

    def save(self):
        """Save cluster plot to file"""
        print "SAVING SCREENGRAB"
        fname = "/home/cat/Pictures/screengrab.png"
        
        if fname:
            fname = str(fname) # convert from QString
            image = self.grabFrameBuffer() # defaults to withAlpha=False, makes no difference
            try:
                image.save(fname)
            except Exception as e:
                QtGui.QMessageBox.critical(
                    self.panel, "Error saving file", str(e),
                    QtGui.QMessageBox.Ok, QtGui.QMessageBox.NoButton)
            print('cluster plot saved to %r' % fname)


    def wheelEvent(self, event):
        self.zoom(event.delta() / 1000)
        self.updateGL()

    '''
    # this specifies a display list, which is sent once, compiled, and then simply referenced
    # later every time the display needs to be updated. However, display lists are static once
    # compiled - none of their attributes can be changed
    def makeObject(self):
        genList = GL.glGenLists(1)
        GL.glNewList(genList, GL.GL_COMPILE)

        GL.glBegin(GL.GL_QUADS)

        x1 = +0.06
        y1 = -0.14
        x2 = +0.14
        y2 = -0.06
        x3 = +0.08
        y3 = +0.00
        x4 = +0.30
        y4 = +0.22

        self.quad(x1, y1, x2, y2, y2, x2, y1, x1)
        self.quad(x3, y3, x4, y4, y4, x4, y3, x3)

        self.extrude(x1, y1, x2, y2)
        self.extrude(x2, y2, y2, x2)
        self.extrude(y2, x2, y1, x1)
        self.extrude(y1, x1, x1, y1)
        self.extrude(x3, y3, x4, y4)
        self.extrude(x4, y4, y4, x4)
        self.extrude(y4, x4, y3, x3)

        Pi = 3.14159265358979323846
        NumSectors = 200

        for i in range(NumSectors):
            angle1 = (i * 2 * Pi) / NumSectors
            x5 = 0.30 * math.sin(angle1)
            y5 = 0.30 * math.cos(angle1)
            x6 = 0.20 * math.sin(angle1)
            y6 = 0.20 * math.cos(angle1)

            angle2 = ((i + 1) * 2 * Pi) / NumSectors
            x7 = 0.20 * math.sin(angle2)
            y7 = 0.20 * math.cos(angle2)
            x8 = 0.30 * math.sin(angle2)
            y8 = 0.30 * math.cos(angle2)

            self.quad(x5, y5, x6, y6, x7, y7, x8, y8)

            self.extrude(x6, y6, x7, y7)
            self.extrude(x8, y8, x5, y5)

        GL.glEnd()
        GL.glEndList()

        return genList

    def quad(self, x1, y1, x2, y2, x3, y3, x4, y4):
        self.qglColor(self.trolltechGreen)

        GL.glVertex3d(x1, y1, -0.05)
        GL.glVertex3d(x2, y2, -0.05)
        GL.glVertex3d(x3, y3, -0.05)
        GL.glVertex3d(x4, y4, -0.05)

        GL.glVertex3d(x4, y4, +0.05)
        GL.glVertex3d(x3, y3, +0.05)
        GL.glVertex3d(x2, y2, +0.05)
        GL.glVertex3d(x1, y1, +0.05)

    def extrude(self, x1, y1, x2, y2):
        self.qglColor(self.trolltechGreen.darker(250 + int(100 * x1)))

        GL.glVertex3d(x1, y1, +0.05)
        GL.glVertex3d(x2, y2, +0.05)
        GL.glVertex3d(x2, y2, -0.05)
        GL.glVertex3d(x1, y1, -0.05)
    '''


class draw_primitives(object):
    
    def __init__(self):
        
        pass    
    
    def soma_shape(self):     

        #f = open('quad_sphere_vertices.csv', 'rt')
        f = open('quad_cube_vertices.csv', 'rt')
        
        vertices_list = list(csv.reader(f))
        f.close()
        vertices = []
        for row in vertices_list:
            #print row[0].split(" ")
            vertices.append([float(row[0].split(" ")[1]), float(row[0].split(" ")[2]), float(row[0].split(" ")[3])])
        #print vertices
        self.vertices = vertices
        
        #f = open('quad_sphere_faces.csv', 'rt')
        f = open('quad_cube_faces.csv', 'rt')
        face_list = list(csv.reader(f))
        f.close()
        triangle_faces=[]
        quad_faces=[]
        for row in face_list:
            #print row[0].split(" ")[1:]
            if (len(row[0].split(" ")[1:]))==4:
                quad_faces.append([float(row[0].split(" ")[1]), float(row[0].split(" ")[2]), float(row[0].split(" ")[3]), float(row[0].split(" ")[4])])
            else:
                triangle_faces.append([float(row[0].split(" ")[1]), float(row[0].split(" ")[2]), float(row[0].split(" ")[3])])
        self.triangle_faces = triangle_faces
        self.quad_faces = quad_faces
        #print self.quad_faces
        
        lowest_vertex_triangle = []
        for j in range(len(self.triangle_faces)):
            lowest_vertex_triangle.append(min(self.vertices[int(self.triangle_faces[j][0])-1][1],
                                self.vertices[int(self.triangle_faces[j][1])-1][1],
                                self.vertices[int(self.triangle_faces[j][2])-1][1]))
        self.lowest_vertex_triangle=lowest_vertex_triangle

        lowest_vertex_quad = []
        for j in range(len(self.quad_faces)):
            lowest_vertex_quad.append(min(self.vertices[int(self.quad_faces[j][0])-1][1],
                    self.vertices[int(self.quad_faces[j][1])-1][1],
                    self.vertices[int(self.quad_faces[j][2])-1][1],
                    self.vertices[int(self.quad_faces[j][3])-1][1]))
        self.lowest_vertex_quad=lowest_vertex_quad


    def spheres(self, coords, size):

        print "Generating spheres"
        #print coords

        #Load sphere segment primitives from .csv file 
        self.soma_shape()

        self.points_quads=[]
        self.points_colors=[]
            
        self.vertices = np.array(self.vertices)
        
        cmap = mpl.cm.get_cmap('jet')
            
        for k in range(len(coords)):
        #for k in range(3):
            soma_xyz=coords[k]
            #print "Coord: ", k, soma_xyz
            #Make spheres at each coordinate point using vertices and quads
            for j in range(len(self.quad_faces)):
                self.points_quads.append([   self.vertices[int(self.quad_faces[j][0])-1]*size+soma_xyz,
                                             self.vertices[int(self.quad_faces[j][1])-1]*size+soma_xyz, 
                                             self.vertices[int(self.quad_faces[j][2])-1]*size+soma_xyz,
                                             self.vertices[int(self.quad_faces[j][3])-1]*size+soma_xyz,
                                             ])

                shade =  (j+20.)/26.#(self.lowest_vertex_quad[j]+2.)/4 
                
                sphere_color = cmap(float(k)/len(coords)) #np.array(color_temp)#/len(coords)*(k*.8)+50
                
                self.points_colors.append([[sphere_color[0]*shade,
                                                sphere_color[1]*shade,
                                                sphere_color[2]*shade
                                                ,255.0]]*4)  #Use 1.0 to get transparency

        #print "quads: ",self.points_quads 
        #print "Pts COLORS: ", self.points_colors
        #print "Length points_quads: ", len(self.points_quads)
        #print "Length points_colors: ", len(self.points_colors)

    def vertices():
        pass

    def simFrame(self, frame_scale):

        #Build frame out of vertexes
        self.points_frame=[]
        frame_scale = 1E3*frame_scale

        vertex_0=[-frame_scale,-frame_scale,frame_scale]
        vertex_1=[frame_scale,-frame_scale,frame_scale]
        vertex_2=[frame_scale,-frame_scale,-frame_scale]
        vertex_3=[-frame_scale,-frame_scale,-frame_scale]
        vertex_4=[-frame_scale,frame_scale,frame_scale]
        vertex_5=[frame_scale,frame_scale,frame_scale]
        vertex_6=[frame_scale,frame_scale,-frame_scale]
        vertex_7=[-frame_scale,frame_scale,-frame_scale]
        
        #********OUTER BOX FRAME
        self.points_frame.append(vertex_0)
        self.points_frame.append(vertex_1)
        self.points_frame.append(vertex_1)
        self.points_frame.append(vertex_2)
        self.points_frame.append(vertex_2)
        self.points_frame.append(vertex_3)
        self.points_frame.append(vertex_3)
        self.points_frame.append(vertex_0)

        self.points_frame.append(vertex_4)
        self.points_frame.append(vertex_5)
        self.points_frame.append(vertex_5)
        self.points_frame.append(vertex_6)
        self.points_frame.append(vertex_6)
        self.points_frame.append(vertex_7)
        self.points_frame.append(vertex_7)
        self.points_frame.append(vertex_4)

        self.points_frame.append(vertex_4)
        self.points_frame.append(vertex_0)
        self.points_frame.append(vertex_5)
        self.points_frame.append(vertex_1)
        self.points_frame.append(vertex_6)
        self.points_frame.append(vertex_2)
        self.points_frame.append(vertex_7)
        self.points_frame.append(vertex_3)

        self.colors_frame = CMAP[[12]*36] # uint8; Need colour for every node, not every vertex; i.e. 2 x no. vertices  


