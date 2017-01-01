#from tsf_ptcs_classes import *
#from sequential_firing import *
#from opengl_pyqt_classes import *
#from sta_utils import *

from mouse import *

from lever_pull_utils import *
import filter
import os
import glob
import multiprocessing as mp

import numpy as np
import matplotlib.pyplot as plt
import struct
import string
import scipy

from matplotlib import cm
import matplotlib.animation as animation
from scipy import ndimage


#Initialize pyqtgraph - needs to be done at beginning - or even better in a main conditional below...
from pyqtgraph.Qt import QtCore, QtGui



#******************************************
#************* DATA FILES *****************
#******************************************
home_dir = '/media/cat/12TB/in_vivo/tim/yuki/'

#mice_names =['AQ2', 'AQ3', 'IA1', 'IA2', 'IA3', 'IJ1', 'IJ2']
mice_names =['AQ5']

  
#******************************************
#********** LOAD MICE  ********************
#******************************************

n_sec = 3                   #number of sec to use for DF/F
window = int(30*n_sec)      #number of frames in window used for analysis
n_pixels = 128              #number of pixels for imaging data.

mice = []
for mouse_name in mice_names:
    if (os.path.exists(home_dir+mouse_name+'/'+mouse_name+'.pkl')==True):
        mouse = Mouse(mouse_name, home_dir, n_sec, window, n_pixels)
        mouse = mouse.load_mouse(mouse_name, home_dir)
        #mouse.save_DFF()
        
        mouse.move_tifs()
        
        mice.append(mouse)
        
    else:
        mouse = Mouse(mouse_name, home_dir, n_sec, window, n_pixels)
        mouse.process()
        mouse.move_tifs()

        mice.append(mouse)
