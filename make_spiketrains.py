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

np.set_printoptions(suppress=True)

#file_name = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-11/2015-12-11-12-allelectinthalamus-iso1/unit_49_channel_16_ptp_091.csv'
#file_name = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-15-allelectrodeinthalamus-is0/unit_59_channel_16_ptp_074.csv'
#file_name = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-2/unit_07_channel_16_ptp_225.csv'
#file_name = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-8-9electrodein-iso0/unit_12_channel_16_ptp_154.csv'
#file_name = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-11/2015-12-11-15-allelectinthalamus-iso0/unit_49_channel_16_ptp_091.csv'
#file_name = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-14-allelectrodeinthalamus-is0/unit_59_channel_16_ptp_074.csv'
file_name = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-11/2015-12-11-17-allelectinthalamus-iso0/unit_49_channel_16_ptp_091.csv'
file_name = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/unit_30_channel_16_ptp_277.csv'

spikes = np.random.random(30000)*3715
spikes = np.sort(spikes)

spikes = []
ctr=0
while ctr < 5500:
    spikes.append(ctr)
    ctr+=0.030

print np.around(spikes, )

np.savetxt(file_name, spikes, fmt='%.5f')
