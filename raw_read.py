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
from lever_pull_utils import fast_mask

from threading import *
import time

import sklearn  
from sklearn import manifold
from tsne import bh_sne
import scipy.ndimage as ndimage

from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import matplotlib.gridspec as gridspec


file_name = '/media/cat/12TB/in_vivo/tim/cat/2016_06_09_test/test_laser_imaging'


data = np.fromfile(file_name, dtype=np.int16)

data = data.reshape((-1, 128, 128))

print data[0]
plt.imshow(data[0], cmap=plt.get_cmap('gray'))
plt.show()
