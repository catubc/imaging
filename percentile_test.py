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

def load_generic_mask(main_dir):
    
    n_pixels=128
    generic_mask_file = main_dir + '/genericmask.txt'
    generic_coords = np.loadtxt(generic_mask_file)
        
    generic_mask_indexes=np.zeros((n_pixels,n_pixels))
    for i in range(len(generic_coords)):
        generic_mask_indexes[generic_coords[i][0]][generic_coords[i][1]] = True

    #Mask centreline also
    for x in range(n_pixels/2-7, n_pixels/2+7, 1):
        for k in range(n_pixels):
            generic_mask_indexes[k][x]=True
    
    return generic_mask_indexes

generic_mask_dir = '/media/cat/12TB/in_vivo/tim/yuki/AQ3/'
main_dir = '/media/cat/12TB/in_vivo/tim/yuki/AQ3/tif_files/AQ3am_Mar3_30Hz/AQ3am_Mar3_30Hz_'
#main_dir = '/media/cat/12TB/in_vivo/tim/yuki/AQ3/tif_files/AQ3pm_Mar7_Day3_30Hz/AQ3pm_Mar7_Day3_30Hz_'
#main_dir = '/media/cat/12TB/in_vivo/tim/yuki/AQ3/tif_files/AQ3pm_Mar9_Day5_30Hz/AQ3pm_Mar9_Day5_30Hz_'
#main_dir = '/media/cat/12TB/in_vivo/tim/yuki/AQ3/tif_files/AQ3am_Apr29_Week8_30Hz/AQ3am_Apr29_Week8_30Hz_'


generic_mask_indexes = load_generic_mask(generic_mask_dir)


block = 6               #pick the # frame block to average together.
for k in range(10):
    data = np.load(main_dir+str(k).zfill(4)+'.npy')
    temp_plot1=data
    norm_block = []
    for k in range(0,len(temp_plot1),block):                          #Average every 6 frame block
        norm_block.append(np.ma.average(temp_plot1[k:k+block], axis=0))
    
    ax=plt.subplot(6,1,1)
    data_1d = np.ma.hstack((temp_plot1))
    data_1d = data_1d.reshape(data_1d.shape[0]*data_1d.shape[1])
    bin_width = 1
    x = np.linspace(np.min(data_1d),np.max(data_1d),200)
    temp_hist = np.histogram(data_1d, bins = x)[0]    
    plt.plot(x[:-1], temp_hist, linewidth=2)
    
    val95 = np.percentile(data_1d, 95)
    plt.plot([val95, val95],[0,np.max(temp_hist)], color='blue')
    val99 = np.percentile(data_1d, 99)
    plt.plot([val99, val99],[0,np.max(temp_hist)], color='red')
    val999 = np.percentile(data_1d, 99.9)
    plt.plot([val999, val999],[0,np.max(temp_hist)], color='green')

    val95_neg = -np.percentile(-data_1d, 95)
    plt.plot([val95_neg, val95_neg],[0,np.max(temp_hist)], color='blue')
    val99_neg = -np.percentile(-data_1d, 99)
    plt.plot([val99_neg, val99_neg],[0,np.max(temp_hist)], color='red')
    val999_neg = -np.percentile(-data_1d, 99.9)
    plt.plot([val999_neg, val999_neg],[0,np.max(temp_hist)], color='green')

    
    ax=plt.subplot(6,1,2)
    norm_block = []
    for k in range(0,len(temp_plot1),block):                          #Average every 6 frame block
        norm_block.append(np.ma.average(temp_plot1[k:k+block], axis=0))
    plt.imshow(np.ma.hstack((norm_block)),vmin=-np.max(np.abs(np.ma.hstack((norm_block)))), vmax=np.max(np.abs(np.ma.hstack((norm_block)))))

    ax=plt.subplot(6,1,3)
    norm_block = []
    for k in range(0,len(temp_plot1),block):                          #Average every 6 frame block
        norm_block.append(np.ma.average(temp_plot1[k:k+block], axis=0))
    stack = np.ma.hstack((norm_block))
    stack = np.ndarray.clip(stack, max=val999)
    #stack [stack>val99]=0.
    plt.imshow(stack,vmin=-np.max(np.abs(stack)), vmax=np.max(np.abs(stack)))


    temp_plot1 = fast_mask(data, generic_mask_indexes)
    norm_block = []
    for k in range(0,len(temp_plot1),block):                          #Average every 6 frame block
        norm_block.append(np.ma.average(temp_plot1[k:k+block], axis=0))
    
    ax=plt.subplot(6,1,4)
    data_1d = np.ma.hstack((temp_plot1))
    data_1d = data_1d.reshape(data_1d.shape[0]*data_1d.shape[1])
    bin_width = 1
    x = np.linspace(np.min(data_1d),np.max(data_1d),200)
    temp_hist = np.histogram(data_1d, bins = x)[0]    
    plt.plot(x[:-1], temp_hist, linewidth=2)
    val95 = np.percentile(data_1d, 95)
    plt.plot([val95, val95],[0,np.max(temp_hist)], color='blue')
    val99 = np.percentile(data_1d, 99)
    plt.plot([val99, val99],[0,np.max(temp_hist)], color='red')
    val999 = np.percentile(data_1d, 99.9)
    plt.plot([val999, val999],[0,np.max(temp_hist)], color='green')

    
    ax=plt.subplot(6,1,5)
    norm_block = []
    for k in range(0,len(temp_plot1),block):                          #Average every 6 frame block
        norm_block.append(np.ma.average(temp_plot1[k:k+block], axis=0))
    plt.imshow(np.ma.hstack((norm_block)),vmin=-np.max(np.abs(np.ma.hstack((norm_block)))), vmax=np.max(np.abs(np.ma.hstack((norm_block)))))

    ax=plt.subplot(6,1,6)
    norm_block = []
    for k in range(0,len(temp_plot1),block):                          #Average every 6 frame block
        norm_block.append(np.ma.average(temp_plot1[k:k+block], axis=0))
    stack = np.ma.hstack((norm_block))
    #stack = np.ndarray.clip(stack, max=val999)
    stack [stack==0]=100.
    plt.imshow(stack,vmin=-np.max(np.abs(stack)), vmax=np.max(np.abs(stack)))







    plt.show()

    
