import sys
sys.path.insert(0, '/home/cat/code/ephys')
from tsf_ptcs_classes import *
from distributions import *
from utils_dongsheng import *
from sta_utils import *

from sklearn import decomposition
from sklearn import datasets

from scipy import stats
from scipy.interpolate import interp1d
import scipy.optimize
import scipy.io
import numpy as np
import time, math
import os.path
import glob #Wildcard searching in dir listing
from pylab import *
import struct, array, csv
import pandas as pd
import matplotlib.mlab as mlab
import itertools
import cPickle as pkl
from scipy import interpolate
import scipy.cluster.hierarchy as sch

from mpl_toolkits.mplot3d import Axes3D

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

colors=['blue','green','violet','lightseagreen','lightsalmon','indianred','lightsalmon','pink','darkolivegreen']

main_dir = '/media/cat/8TB/in_vivo/tim/dongsheng/'

file_dirs = []
file_names = []

##########**************** LOAD ALL FILES DATA IN EXPERIMENT DIRECTORY ****************
#file_dirs.append('/media/cat/12TB/in_vivo/tim/dongsheng/2015-5-6/')  #*************************************
#file_names.append([
#'2015-5-6-1',       #cortex
#'2015-5-6-2',       #cortex
#'2015-5-6-3',      #subcortex
#'2015-5-6-4',      #subcortex
#'2015-5-6-5',      #subcortex
#'2015-5-6-6'       #subcortex

##'2015-5-6-7-FLstim',
##'2015-5-6-8-hlstim',
##'2015-5-6-9-v1',
##'2015-5-6-10-A1',
#])

#file_dirs.append('/media/cat/12TB/in_vivo/tim/dongsheng/2015-7-22/')  #*************************************
#file_names.append([
#'2015-7-22-1',     #cortex
#'2015-7-22-2',     #cortex
#'2015-7-22-3',     #cortex
#'2015-7-22-4',     #subcortex
#'2015-7-22-5',     #subcortex
#'2015-7-22-6',     #subcortex
#'2015-7-22-7',     #subcortex
#'2015-7-22-8'      #subcortex

##'2015-7-22-9-W1', 
##'2015-7-22-10-w2', 
##'2015-7-22-11-w3',
##'2015-7-22-12-v1',
##'2015-7-22-13-V2',
##'2015-7-22-14-FL',
##'2015-7-22-15-FL',
##'2015-7-22-16-HL',
##'2015-7-22-17-HL',
##'2015-7-22-18-A1',
##'2015-7-22-19-A2'
#])

#file_dirs.append('/media/cat/12TB/in_vivo/tim/dongsheng/2015-7-23/')  #*************************************
#file_names.append([
#'2015-7-23-1',     #cortex
#'2015-7-23-2',     #subcortex
#'2015-7-23-3',     #subcortex  
#'2015-7-23-16'     #subcortex

##'2015-7-23-5-w1',
##'2015-7-23-6-w2',
##'2015-7-23-7-v1',
##'2015-7-23-8-v2',
##'2015-7-23-9-a1',
##'2015-7-23-10-A2',
##'2015-7-23-11-FL',
##'2015-7-23-12-FL',
##'2015-7-23-13-HL',
##'2015-7-23-14-HL',
#])

#file_dirs.append('/media/cat/12TB/in_vivo/tim/dongsheng/2015-7-30/')  #*************************************
#file_names.append([
#'2015-7-30-4',
#'2015-7-30-5', 
#'2015-7-30-6', 
#'2015-7-30-7'
#])

file_dirs.append('/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/') #*************************************
file_names.append([
'2015-11-18-1-9electrodein', 
'2015-11-18-2-9electrodein', 
'2015-11-18-3-9electrodein',
'2015-11-18-4-9electrodein-iso0.5',
'2015-11-18-5-9electrodein-iso0.5',
'2015-11-18-6-9electrodein-iso0.5',
'2015-11-18-7-9electrodein-iso0',
'2015-11-18-8-9electrodein-iso0',
'2015-11-18-9-9electrodein-iso0',
'2015-11-18-10-allelectrodein-iso0',
'2015-11-18-11-allelectrodein-iso0',
'2015-11-18-12-allelectrodein-iso0',
'2015-11-18-13-deep-iso0',
'2015-11-18-14-deep-iso0', 
'2015-11-18-15-deep-iso0', 
'2015-11-18-16-deep-iso0', 
'2015-11-18-17-deep-iso0', 
'2015-11-18-18-deep-iso0'
])


file_dirs.append('/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/') #*************************************
file_names.append([
'2015-11-27-1-10electrodein-iso1.0',
'2015-11-27-2-10electrodein-iso1', 
'2015-11-27-4-10electrodein-iso0', 
'2015-11-27-5-10electrodein-iso0',
'2015-11-27-7-16electrodein-iso1',
'2015-11-27-8-16electrodein-iso1',
'2015-11-27-10-16electrodein-iso0',
'2015-11-27-11-16electrodein-iso0',
'2015-11-27-13-deep-iso1',
'2015-11-27-14-deep-iso1',
'2015-11-27-16-deep-iso0',

##'2015-11-27-2-10electrodein-iso1_mua', 
##'2015-11-27-5-10electrodein-iso0_mua',
##'2015-11-27-10-16electrodein-iso0_mua',
##'2015-11-27-14-deep-iso1_mua',
##'2015-11-27-16-deep-iso0_mua'

##'2015-11-27-2-10electrodein-iso1_lfp', 
##'2015-11-27-5-10electrodein-iso0_lfp',
##'2015-11-27-14-deep-iso1_lfp',
##'2015-11-27-16-deep-iso0_lfp'
])

file_dirs.append('/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/') #*************************************
file_names.append([
'2015-12-1-1-10electrodeiniso1',
'2015-12-1-2-10electrodeiniso1',
'2015-12-1-3-10electrodeiniso1',
'2015-12-1-5-10electrodeiniso0',
'2015-12-1-6-10electrodeiniso0',
'2015-12-1-8-allelectrodeiniso0.8', 
'2015-12-1-9-allelectrodeiniso0.8', 
'2015-12-1-11-allelectrodeiniso0', 
'2015-12-1-12-allelectrodeiniso0',
'2015-12-1-14-5electrodeinthalamus-iso0.8',
'2015-12-1-15-5electrodeinthalamus-iso0.8',
'2015-12-1-18-5electrodeinthalamus-iso0',   #SINGLE
'2015-12-1-19-allelectrodeinthalamus-iso0.8',  #SINGLE
'2015-12-1-20-allelectrodeinthalamus-iso0.8',
'2015-12-1-22-allelectrodeinthalamus-iso0',
'2015-12-1-23-allelectrodeinthalamus-iso0', 
'2015-12-1-24-allelectrodeinthalamus-iso0'

#'2015-12-1-3-10electrodeiniso1_mua',
#'2015-12-1-6-10electrodeiniso0_mua',
#'2015-12-1-20-allelectrodeinthalamus-iso0.8_mua',
#'2015-12-1-23-allelectrodeinthalamus-iso0_mua', 

#'2015-12-1-3-10electrodeiniso1_lfp',
#'2015-12-1-6-10electrodeiniso0_lfp',
#'2015-12-1-20-allelectrodeinthalamus-iso0.8_lfp',
#'2015-12-1-23-allelectrodeinthalamus-iso0_lfp', 
])

file_dirs.append('/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/') #*************************************
file_names.append([
#'2015-12-2-2-10electrodeincortex-iso1', 
'2015-12-2-3-10electrodeincortex-iso1',
'2015-12-2-4-10electrodeincortex-iso1',
'2015-12-2-6-10electrodeincortex-iso0',
'2015-12-2-7-10electrodeincortex-iso0',
'2015-12-2-9-allelectrodeincortex-iso1',  #SINGLE
'2015-12-2-11-allelectrodeinthalamus-iso1',
'2015-12-2-12-allelectrodeinthalamus-iso1',
'2015-12-2-14-allelectrodeinthalamus-is0',
'2015-12-2-15-allelectrodeinthalamus-is0',

#'2015-12-2-3-10electrodeincortex-iso1_mua',
#'2015-12-2-6-10electrodeincortex-iso0_mua',
#'2015-12-2-12-allelectrodeinthalamus-iso1_mua',
#'2015-12-2-14-allelectrodeinthalamus-is0_mua',

#'2015-12-2-3-10electrodeincortex-iso1_lfp',
#'2015-12-2-6-10electrodeincortex-iso0_lfp',
#'2015-12-2-12-allelectrodeinthalamus-iso1_lfp',
#'2015-12-2-14-allelectrodeinthalamus-is0_lfp'
])

file_dirs.append('/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-11/') #*************************************
file_names.append([
'2015-12-11-1-10electincortex-iso1.0',
'2015-12-11-2-10electincortex-iso1',
'2015-12-11-3-10electincortex-iso1',
'2015-12-11-5-10electincortex-iso0',
'2015-12-11-7-10electincortex-iso1',
'2015-12-11-8-10electincortex-iso1',
'2015-12-11-9-allelectincortex',
'2015-12-11-10-5electinthalamusorcpu',
'2015-12-11-11-allelectinthalamus-iso1',
'2015-12-11-12-allelectinthalamus-iso1',
'2015-12-11-13-allelectinthalamus-iso1',
'2015-12-11-15-allelectinthalamus-iso0',
'2015-12-11-16-allelectinthalamus-iso0',
'2015-12-11-17-allelectinthalamus-iso0',
'2015-12-11-18-allelectinthalamus-iso1'

#'2015-12-11-3-10electincortex-iso1_mua',
#'2015-12-11-5-10electincortex-iso0_mua',
#'2015-12-11-12-allelectinthalamus-iso1_mua',
#'2015-12-11-16-allelectinthalamus-iso0_mua',

#'2015-12-11-3-10electincortex-iso1_lfp',
#'2015-12-11-5-10electincortex-iso0_lfp',
#'2015-12-11-12-allelectinthalamus-iso1_lfp',
#'2015-12-11-16-allelectinthalamus-iso0_lfp',
])

file_dirs.append('/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-16/') #*************************************
file_names.append([
#'2015-12-16-1-11electrodeincortex-iso1',
'2015-12-16-2-11electrodeincortex-iso1',
'2015-12-16-4-11electrodeincortex-iso0',
'2015-12-16-5-11electrodeincortex-iso0',
'2015-12-16-7-11electrodeincortex-iso1',
'2015-12-16-8-4electrodeincpu-iso1',
'2015-12-16-9-allelectrodeinZI-iso1',
'2015-12-16-10-allelectrodeinZI-iso1',
'2015-12-16-11-allelectrodeinZI-iso1',
'2015-12-16-13-allelectrodeinZI-iso0',
'2015-12-16-14-allelectrodeinZI-iso0',
'2015-12-16-16-allelectrodeinZI-iso1',
'2015-12-16-17-allelectrodeinZI-iso1'
])


file_dirs.append('/media/cat/8TB/in_vivo/tim/dongsheng/2016-1-11/') #*************************************
file_names.append([
'2016-1-11-1-10electrodeincortex-iso1.2',
'2016-1-11-2-10electrodeincortex-iso1',
'2016-1-11-3-10electrodeincortex-iso1',
'2016-1-11-5-10electrodeincortex-iso0',
'2016-1-11-6-10electrodeincortex-iso0',
'2016-1-11-8-10electrodeincortex-iso1',
'2016-1-11-9-allelectrodeincortex-iso1',
'2016-1-11-10-13electrodeinthalamus-iso1',
'2016-1-11-11-13electrodeinthalamus-iso1',
'2016-1-11-12-13electrodeinthalamus-iso1',
'2016-1-11-14-13electrodeinthalamus-iso0',
'2016-1-11-15-13electrodeinthalamus-iso0',
'2016-1-11-17-13electrodeinthalamus-iso1',
'2016-1-11-18-13electrodeinthalamus-iso1'
])

file_dirs.append('/media/cat/8TB/in_vivo/tim/dongsheng/2016-1-14/') #*************************************
file_names.append([
'2016-1-14-1-10electrodein cortex-iso1',
'2016-1-14-2-10electrodein cortex-iso1',
'2016-1-14-3-10electrodein cortex-iso1',
'2016-1-14-5-10electrodein cortex-iso0',
'2016-1-14-6-10electrodein cortex-iso0',
'2016-1-14-7-10electrodein cortex-iso0',
'2016-1-14-9-10electrodein cortex-iso1',
'2016-1-14-10-allelectrodein cortex-iso1',
'2016-1-14-11-allelectrodeinthalamus-iso1',
'2016-1-14-12-allelectrodeinthalamus-iso1',
'2016-1-14-13-allelectrodeinthalamus-iso1',
'2016-1-14-15-allelectrodeinthalamus-iso0',
'2016-1-14-16-allelectrodeinthalamus-iso0',
'2016-1-14-17-allelectrodeinthalamus-iso0'
])


#*********************** PARAMETERS *******************
area_names = ['hindlimb', 'forelimb', 'barrel','retrosplenial','visual', 'motor', 'pta', 'acc'] 

#state_match = "anesthetized"
state_match = "awake"

depth_match = "cortex"
#depth_match = "subcortical"

data = Object_temp()

#****************** PARCELATION CLUSTERING ************
if False: parcelation_clustering()

#****************** ROI CLUSTERING ********************
if True: roi_clustering(data, file_dirs, file_names, main_dir, state_match, depth_match)

print len(data.cortex)
print len(data.subcortical)


#**************** PLOT STACKED VECTOR DATA ************
if True:
    
    #fig, ax = plt.subplots(nrows=1, ncols=2)
    fig = plt.figure()
    gs = gridspec.GridSpec(12,8)
    #ax = plt.subplot(gs[0,0:1])

    #Reorder each ROI area by location of activation
    data1 = data.cortex
    if len(data1)>0:
        max_val = []
        for j in range(len(data1)):  #Loop over all cells
            max_val.append(np.argmax(np.abs(data1[j][466:600])))
            #max_val.append(np.argmax(data1[j][466:600]))
            
        indexes = np.argsort(max_val)
        img_out1= []
        counter=0
        reassigned_cell = []
        for p in indexes:
            if (state_match=='awake') and ((counter>210) or (counter<5)):# and (data.cortex_ch[p]>10):  #Reassign cortical cells to subcortex
                reassigned_cell.append(data1[p])
            else: 
                if (state_match=='anesthetized') and (counter>250) and (data.cortex_ch[p]>10):  #Reassign cortical cells to subcortex
                    reassigned_cell.append(data1[p])
                else:
                    if (np.max(data1[p])>-np.min(data1[p])):
                        temp = data1[p]/np.max(data1[p])
                    else:
                        temp = data1[p]/-np.min(data1[p])
                        
                    img_out1.append(temp)
                    
            if counter>200: print data.cortex_ch[p]
            
            counter+=1

        #data.subcortical = np.concatenate((data.subcortical, reassigned_cell), axis=0)

        img_out1 = np.array(img_out1)
        v_min = np.min(img_out1); v_max = np.max(img_out1)
        if v_max<(-v_min): v_max=-v_min
        else: v_min = -v_max

        for p in range(8):
            ax1 = fig.add_subplot(gs[0:4, p])

            #plt.imshow(img_out1[:,p*200:(p+1)*200],  vmin = v_min, vmax = v_max, cmap=plt.get_cmap('RdBu_r'), aspect='auto', interpolation='none')
            plt.imshow(img_out1[:,p*200:(p+1)*200],  vmin = v_min, vmax = v_max, aspect='auto', interpolation='none')
            
            xx = np.arange(100,100+len(area_names)*200,200)
            ax1.xaxis.tick_top()
            for k in range(len(area_names)-1):
                plt.plot([200+k*200,200+k*200],[-0.5,len(img_out1)], linewidth=3, color='black')
            for j in range(len(area_names)):
                plt.plot([100.+200*j,100.+200*j],[-0.5,len(img_out1)], 'r--', linewidth = 2, color='black', alpha=0.5)

            plt.ylim(len(img_out1)-0.5,-0.5)
            plt.title(area_names[p])
            plt.xlim(0,200)
            if p>0: ax1.get_yaxis().set_visible(False)
            ax1.get_xaxis().set_visible(False)

    #Plot histograms
    for k in range(8):
        ax1 = fig.add_subplot(gs[4, k])
        img_avg = np.average(img_out1[:,k*200:(k+1)*200],axis=0)
        temp_std = np.std(img_out1[:,k*200:(k+1)*200], axis=0)
        xx = np.linspace(0,200,200)
        ax1.fill_between(xx, img_avg+temp_std, img_avg-temp_std, facecolor='black', alpha=0.3)
                
        plt.plot(xx,img_avg, color='black', linewidth=3)
        plt.ylim(-0.5,1.0)
        plt.plot([0,200], [0.,0.], 'r--', color='black')
        plt.plot([100,100], [-0.5,1.], 'r--', color='black')
        ax1.get_yaxis().set_visible(False)
        ax1.get_xaxis().set_visible(False)

    
    #************** Plot Subcortical Data ********
    #Reorder each ROI area by location of activation
    data1 = data.subcortical
    if len(data1)>0:
        max_val = []
        for j in range(len(data1)):  #Loop over all cells
            max_val.append(np.argmax(np.abs(data1[j][466:600])))
            #max_val.append(np.argmax(data1[j][466:600]))

        indexes = np.argsort(max_val)
        img_out2= []
        for p in indexes:
            if (np.max(data1[p])>-np.min(data1[p])):
                temp = data1[p]/np.max(data1[p])
            else:
                temp = data1[p]/-np.min(data1[p])
                
            img_out2.append(temp)


        img_out2 = np.array(img_out2)
        v_min = np.min(img_out2); v_max = np.max(img_out2)
        if v_max<(-v_min): v_max=-v_min
        else: v_min = -v_max

        for p in range(8):
            ax1 = fig.add_subplot(gs[6:10, p])

            #im = ax1.imshow(img_out2[:,p*200:(p+1)*200], vmin=v_min, vmax=v_max, cmap=plt.get_cmap('RdBu_r'), aspect='auto', interpolation='none')
            im = ax1.imshow(img_out2[:,p*200:(p+1)*200], vmin=v_min, vmax=v_max, aspect='auto', interpolation='none')
            xx = np.arange(100,100+len(area_names)*200,200)
            ax1.xaxis.tick_top()
            for k in range(len(area_names)-1):
                plt.plot([200+k*200,200+k*200],[-0.5,len(img_out2)], linewidth=3, color='black')
            for j in range(len(area_names)):
                plt.plot([100.+200*j,100.+200*j],[-0.5,len(img_out2)], 'r--', linewidth = 2, color='black', alpha=0.5)
            #plt.ylabel("Cell id (max/min order)")
            #plt.xticks(xx, area_names) 
            plt.ylim(len(img_out2)-0.5,-0.5)
            plt.title(area_names[p])
            #plt.xlim(0,len(img_out1[0]))
            plt.xlim(0,200)
            if p>0: ax1.get_yaxis().set_visible(False)
            ax1.get_xaxis().set_visible(False)
            
        plt.suptitle("All cells: "+ state_match)
        
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.75, 0.05, 0.1, 0.105])
        cbar = fig.colorbar(im,orientation='horizontal', ticks=[-1, 0, 1], cax=cbar_ax)
        cbar.ax.set_xticklabels(['-1', '0', '1'], fontsize=15)  # horizontal colorbar

    #Plot histograms
    for k in range(8):
        ax1 = fig.add_subplot(gs[10, k])
        img_avg = np.average(img_out2[:,k*200:(k+1)*200],axis=0)

        temp_std = np.std(img_out2[:,k*200:(k+1)*200], axis=0)
        xx = np.linspace(0,200,200)
        ax1.fill_between(xx, img_avg+temp_std, img_avg-temp_std, facecolor='black', alpha=0.3)

        plt.plot(xx, img_avg, color='black', linewidth=3)
        plt.ylim(-0.5,1.0)
        plt.plot([0,200], [0.,0.], 'r--', color='black')
        plt.plot([100,100], [-0.5,1.], 'r--', color='black')
        ax1.get_yaxis().set_visible(False)
        ax1.get_xaxis().set_visible(False)
        
    plt.show()

#********* Plot correlation matrices for cells:
if True: 
    data_in = data.subcortical
    img_ave = []
    for i in range(len(data_in)):
        ax = plt.subplot(6,6,i+1)
        
        #*******Plot correlation matrix
        temp_array = np.split(data_in[i],8)
        img_corr = np.corrcoef(temp_array)
        #plt.imshow(img_corr, origin='lower', interpolation='none', cmap=plt.get_cmap('jet'))
        img_ave.append(img_corr)

        #***** Plot clustered corr matrix
        D_cortex = np.corrcoef(np.split(data_in[i],8))
        Y = sch.linkage(D_cortex, method='centroid')
        Z1_cortex = sch.dendrogram(Y, orientation='bottom')
        Y = sch.linkage(D_cortex, method='single')
        Z2_cortex = sch.dendrogram(Y)

        #Cortex plot
        idx1 = Z1_cortex['leaves']
        idx2 = Z2_cortex['leaves']
        D = D_cortex[idx1,:]
        D = D[:,idx2]
        im = ax.matshow(D, origin='lower',interpolation='none', cmap=plt.get_cmap('jet'))
        #plt.title(depth_match + " " + state_match + "\n(in cluster order)")
        #axmatrix.set_xticks(idx2)
        #axmatrix.set_yticks(idx1)

        xx = np.linspace(0,len(D_cortex)-1,len(D_cortex))
        label_idx2 = []
        for k in range(len(idx2)):
            label_idx2.append(str(area_names[idx2[k]][0:1]))
        label_idx1 = []
        for k in range(len(idx1)):
            label_idx1.append(str(area_names[idx1[k]][0:1]))        
        
        plt.xticks(xx, label_idx2)
        plt.xticks(rotation=90)
        plt.yticks(xx, label_idx1)

        plt.ylim(-.5,-.5+len(D_cortex))
        plt.xlim(-.5,-.5+len(D_cortex))

    
    #***** Plot average total matrix
    ax = plt.subplot(6,6,31)
    img_avg = np.average(np.array(img_ave), axis=0)
    
    D_cortex = np.corrcoef(img_avg)
    Y = sch.linkage(D_cortex, method='centroid')
    Z1_cortex = sch.dendrogram(Y, orientation='bottom')
    Y = sch.linkage(D_cortex, method='single')
    Z2_cortex = sch.dendrogram(Y)

    #Cortex plot
    idx1 = Z1_cortex['leaves']
    idx2 = Z2_cortex['leaves']
    D = D_cortex[idx1,:]
    D = D[:,idx2]
    im = ax.imshow(D, origin='lower',interpolation='none', cmap=plt.get_cmap('jet'))
    #plt.title(depth_match + " " + state_match + "\n(in cluster order)")
    #axmatrix.set_xticks(idx2)
    #axmatrix.set_yticks(idx1)

    xx = np.linspace(0,len(D_cortex)-1,len(D_cortex))
    label_idx2 = []
    for k in range(len(idx2)):
        label_idx2.append(str(area_names[idx2[k]][0:1]))
    label_idx1 = []
    for k in range(len(idx1)):
        label_idx1.append(str(area_names[idx1[k]][0:1]))        
    
    plt.xticks(xx, label_idx2)
    plt.xticks(rotation=90)
    plt.yticks(xx, label_idx1)

    plt.ylim(-.5,-.5+len(D_cortex))
    plt.xlim(-.5,-.5+len(D_cortex))

    plt.show()

#convert matrices to vectors
if False: vectors = parcelation_1D(main_dir)

#***************** PLOT DENDROGRAMS *******************
if False: plot_dendrograms()

#***************** CLUSTER DATA ****************
if False: 
    
    #from sequential_firing import multi_dim_scaling
    #data = multi_dim_scaling(vectors,2)

    labels = KMEANS(vectors,8)

    print labels.shape

    #Remove smallest clusters
    labels_unique = np.unique(labels, return_counts=True)[1]
    print labels_unique

    del_elements = []
    for i in range(len(labels_unique)):
        if labels_unique[i]<5:
            del_elements.extend(np.where(labels==i)[0])

    del_elements = np.sort(del_elements)[::-1]
    for i in range(len(del_elements)):
        data = scipy.delete(data, del_elements[i], 0)

    print data.shape

    labels = KMEANS(data,8)

    #Remove smallest clusters
    labels_unique = np.unique(labels, return_counts=True)[1]
    print labels_unique

    del_elements = []
    for i in range(len(labels_unique)):
        if labels_unique[i]<5:
            del_elements.extend(np.where(labels==i)[0])

    del_elements = np.sort(del_elements)[::-1]
    for i in range(len(del_elements)):
        data = scipy.delete(data, del_elements[i], 0)

    print data.shape

    labels = KMEANS(data,6)
        
    
#****************** MEAN SHIFT **************
if False:  labels = MEANSHIFT()

#**************** PLOT 3D DATA **************
if False: plot_3D(labels)




quit()

#PLot cortex ROI
ax1 = plt.subplot(221)
bins = np.linspace(0,len(area_names),len(area_names)+1)
yy =  np.histogram(left_tag_cortex_excitatory, bins = bins)[0]
max_yy  = max(yy)
for i in range(len(area_names)):
    plt.bar(i, yy[i], 1.0, color=colors[i], alpha=0.65)
original_xx = np.arange(0.5,0.5+len(area_names),1)
xx = area_names
plt.xticks(original_xx, xx)
plt.title("Cortex - Excitatory Connections")

#Plot subcortical ROI dynamics #NEED TO GET MAX HISTOGRAM VALUE FIRST
bins = np.linspace(0,len(area_names),len(area_names)+1)
yy =  np.histogram(left_tag_cortex_inhibitory, bins = bins)[0]
max_yy = max(max_yy, max(yy))
plt.ylim(0,max_yy)


ax1 = plt.subplot(222)
for i in range(len(area_names)):
    plt.bar(i, yy[i], 1.0, color=colors[i], alpha=0.65)
original_xx = np.arange(0.5,0.5+len(area_names),1)
xx = area_names
plt.ylim(0,max_yy)
plt.xticks(original_xx, xx)
plt.title("Cortex - Inhibitory Connections")

#PLot cortex ROI
ax1 = plt.subplot(223)
bins = np.linspace(0,len(area_names),len(area_names)+1)
yy =  np.histogram(left_tag_subcortical_excitatory, bins = bins)[0]
max_yy  = max(yy)
for i in range(len(area_names)):
    plt.bar(i, yy[i], 1.0, color=colors[i], alpha=0.65)
original_xx = np.arange(0.5,0.5+len(area_names),1)
xx = area_names
plt.xticks(original_xx, xx)
plt.title("Subcortical - Excitatory Connections")

#Plot subcortical ROI dynamics
bins = np.linspace(0,len(area_names),len(area_names)+1)
yy =  np.histogram(left_tag_subcortical_inhibitory, bins = bins)[0]
max_yy = max(max_yy, max(yy))
plt.ylim(0,max_yy)

ax1 = plt.subplot(224)
for i in range(len(area_names)):
    plt.bar(i, yy[i], 1.0, color=colors[i], alpha=0.65)
original_xx = np.arange(0.5,0.5+len(area_names),1)
xx = area_names
plt.xticks(original_xx, xx)
plt.title("Subcortical - Inhibitory Connections")
plt.ylim(0,max_yy)



plt.suptitle("All mice location of highest network relationship \n"+"Anesthetic state: "+state_match, fontsize=20)
plt.show()
                
                               
                
                







##PLOT Distance matrix

## Compute and plot first dendrogram.
#D_cortex = np.corrcoef(total_cortex_all)
#Y = sch.linkage(D_cortex, method='centroid')
#Z1_cortex = sch.dendrogram(Y, orientation='bottom')
#Y = sch.linkage(D_cortex, method='single')
#Z2_cortex = sch.dendrogram(Y)

#D_subcortical = np.corrcoef(total_subcortical_all)
#Y = sch.linkage(D_subcortical, method='centroid')
#Z1_subcortical = sch.dendrogram(Y, orientation='bottom')
#Y = sch.linkage(D_subcortical, method='single')
#Z2_subcortical = sch.dendrogram(Y)
#plt.close()


##Cortex plot
#ax1 = plt.subplot(241)
#img_out = np.corrcoef(total_cortex_all)
#plt.imshow(img_out, origin='lower', interpolation='none', cmap=plt.get_cmap('jet'))
#n_cells = 0
#for i in range(len(mouse_cortex_counter)):
    #n_cells += mouse_cortex_counter[i]
    #plt.plot([-0.5,len(total_cortex_all)-0.5], [n_cells-0.5,n_cells-0.5],'r--', color='black', linewidth=3)
#plt.ylim(0,len(total_cortex_all))
#plt.xlim(0,len(total_cortex_all))
#plt.title("Cortex " + state_match + "\n (in order of experiments)")

##Distance clusters
#axmatrix = plt.subplot(242)
#idx1 = Z1_cortex['leaves']
#idx2 = Z2_cortex['leaves']
#D = D_cortex[idx1,:]
#D = D[:,idx2]
#im = axmatrix.matshow(D, origin='lower',interpolation='none', cmap=plt.get_cmap('jet'))
#plt.title("Cortex " + state_match + "\n(in cluster order)")
#plt.ylim(0,len(total_cortex_all))
#plt.xlim(0,len(total_cortex_all))
#axmatrix.set_xticks([])
#axmatrix.set_yticks([])

#axmatrix = plt.subplot(243)
#im = axmatrix.matshow(D, origin='lower',vmin=0, interpolation='none', cmap=plt.get_cmap('jet'))
#plt.title("Cortex " + state_match + "\n(high correlations)")
#plt.ylim(0,len(total_cortex_all))
#plt.xlim(0,len(total_cortex_all))
#axmatrix.set_xticks([])
#axmatrix.set_yticks([])

#axmatrix = plt.subplot(244)
#im = axmatrix.matshow(-D, origin='lower',vmin=0, interpolation='none', cmap=plt.get_cmap('jet'))
#plt.title("Cortex " + state_match + "\n(high anticorrelations)")
#plt.ylim(0,len(total_cortex_all))
#plt.xlim(0,len(total_cortex_all))
#axmatrix.set_xticks([])
#axmatrix.set_yticks([])


##********************************Subcortical plot
#ax1 = plt.subplot(245)
#img_out = np.corrcoef(total_subcortical_all)
#plt.imshow(img_out, origin='lower', interpolation='none',cmap=plt.get_cmap('jet'))
#n_cells = 0
#for i in range(len(mouse_subcortical_counter)):
    #n_cells += mouse_subcortical_counter[i]
    #plt.plot([-0.5,len(total_subcortical_all)-0.5], [n_cells-0.5,n_cells-0.5],'r--', color='black', linewidth=3)
#plt.title("Subcortical " + state_match + "\n (in order of experiments)")
#plt.ylim(0,len(total_subcortical_all))
#plt.xlim(0,len(total_subcortical_all))


##Distance clusters
#axmatrix = plt.subplot(246)
#idx1 = Z1_subcortical['leaves']
#idx2 = Z2_subcortical['leaves']
#print idx1
#print idx2

#D = D_subcortical[idx1,:]
#D = D[:,idx2]

#im = axmatrix.matshow(D, origin='lower', interpolation='none',cmap=plt.get_cmap('jet'))
#plt.title("Subcortical " + state_match + "\n(in cluster order)")
#plt.ylim(0,len(total_subcortical_all))
#plt.xlim(0,len(total_subcortical_all))
#axmatrix.set_xticks([])
#axmatrix.set_yticks([])


#axmatrix = plt.subplot(247)
#im = axmatrix.matshow(D, origin='lower', vmin=0, interpolation='none',cmap=plt.get_cmap('jet'))
#plt.title("Subcortical " + state_match + "\n(high correlations)")
#plt.ylim(0,len(total_subcortical_all))
#plt.xlim(0,len(total_subcortical_all))
#axmatrix.set_xticks([])
#axmatrix.set_yticks([])


#axmatrix = plt.subplot(248)
#im = axmatrix.matshow(-D, origin='lower', vmin=0, interpolation='none',cmap=plt.get_cmap('jet'))
#plt.title("Subcortical " + state_match + "\n(high anticorrelations)")
#plt.ylim(0,len(total_subcortical_all))
#plt.xlim(0,len(total_subcortical_all))
#axmatrix.set_xticks([])
#axmatrix.set_yticks([])

#plt.show()

#axmatrix = plt.subplot(111)
#idx1 = Z1_subcortical['leaves']
#idx2 = Z2_subcortical['leaves']
#print idx1
#print idx2

#D = D_subcortical[idx1,:]
#D = D[:,idx2]

#im = axmatrix.matshow(D, origin='lower', interpolation='none',cmap=plt.get_cmap('jet'))
#plt.title("Subcortical " + state_match + "\n(in cluster order)")
#plt.ylim(0,len(total_subcortical_all))
#plt.xlim(0,len(total_subcortical_all))

#xx = np.linspace(0,len(D),len(D))
#plt.xticks(xx, idx2)
#plt.xticks(rotation=90)

#plt.yticks(xx, idx1)

##axmatrix.set_yticks(idx2)
##plt.yticks(rotation=90)
#axmatrix.tick_params(axis='both', which='both', labelsize=8)

#plt.show()





