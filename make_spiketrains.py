import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.ndimage as ndimage
import scipy
import os
from mpl_toolkits.mplot3d import Axes3D
from sklearn import *
import glob
import math
from math import pi
from sklearn import decomposition
from sklearn import datasets
from scipy import interpolate
import scipy.cluster.hierarchy as sch
from scipy.signal import butter, filtfilt, cheby1
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg
import pyqtgraph.opengl as gl

file_name = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-11/2015-12-11-12-allelectinthalamus-iso1/unit_49_channel_16_ptp_091.csv'
np.set_printoptions(suppress=True)

spikes = np.random.random(30000)*3710
spikes = np.sort(spikes)
print np.around(spikes, )

np.savetxt(file_name, spikes, fmt='%.5f')
