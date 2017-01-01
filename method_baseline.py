#Check different baseline removal methods Dongsheng/Tim GCamp data

from tsf_ptcs_classes import *
from distributions import *
from sequential_firing import *

from scipy import stats
from scipy import signal

import numpy as np
import time, math
import sys
import os.path
import multiprocessing as mp
import glob #Wildcard searching in dir listing

from pylab import *
from scipy.interpolate import interp1d
import struct, array, csv
import scipy.optimize 
from scipy.optimize import curve_fit
import pandas as pd
from scipy.interpolate import UnivariateSpline
import matplotlib.mlab as mlab
from matplotlib import animation
from PIL import Image

from sta_utils import *

file_dirs = []
file_names = []
#########**************** LOAD ALL FILES DATA IN EXPERIMENT DIRECTORY ****************
#file_dirs.append('/media/cat/12TB/in_vivo/tim/2015-5-6/')  #*************************************
#file_names.append([
#'2015-5-6-2',       #cortex
#])


file_dirs.append('/media/cat/12TB/in_vivo/tim/2015-12-2/') #*************************************
file_names.append([
#'2015-12-2-3-10electrodeincortex-iso1',
'2015-12-2-7-10electrodeincortex-iso0',
])

colors = ['red', 'blue', 'green']
plot_strings = ['all', 'burst','1sec']

window = 3

dir_counter = 0
for file_dir in file_dirs:
    print file_dir
    print file_names[dir_counter]

    for file_name in file_names[dir_counter]:
        #**************** Manual override temp_names for single run *********************
        files = os.listdir(file_dir+file_name)
        temp_names = []
        for file_ in files:
            if ("unit_" in file_):
                if ('map' in file_):
                    pass
                else:
                    temp_names.append(file_)
        
        #Convert unit spike file names to integers
        units = []
        channels = []
        ptps = []
        for file_ in temp_names:
            #print files
            units.append(int(file_[5:7]))
            channels.append(int(file_[16:18]))
            ptps.append(int(file_[23:26]))
            
        #fig = plt.
        for unit in units:
            counter = 1
            img=[]
            for plot_string in plot_strings:
                
                #Plot static maps
                map_name = glob.glob(file_dir + file_name+'/'+file_name+'_'+plot_string+'_maxmap_unit_'+str(unit).zfill(2)+'*.npy')
                img.append(np.load(map_name[0]))
                ax=plt.subplot(3,3,counter)
                plt.imshow(img[counter-1])
                n_spikes = map_name[0][-9:-4]
                plt.title(plot_string+ "\n#spks: "+str(n_spikes), fontsize=10)
                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)
                if counter==1:
                    ax.get_yaxis().set_visible(True)
                    plt.ylabel("Normalized DF/F", fontsize=10)

                
                #Plot time courses
                tc_name = glob.glob(file_dir + file_name+'/time_course_data_'+file_name+'_'+plot_string+'_'+str(window)+'sec_window_unit'+str(unit).zfill(2)+'*')
                traces=[]
                with open(tc_name[0], "r") as f:
                    data = csv.reader(f)
                    data.next()
                    data.next()
                    for row in data:
                        traces.append([float(i) for i in row])
               
                #print len(traces)
               
                ax=plt.subplot(3,3,6+counter)
                frame_width = float(len(traces[29]))/(window*2.)
                xx = np.arange(-window,window,1/frame_width)
                plt.plot(xx, traces[29],color=colors[counter-1], linewidth=3)
                plt.ylim(-0.02,0.05)
                plt.xlabel("Time (sec)", fontsize=10)
                plt.plot([0,0],[-1,1], color='black')
                plt.plot([-3,3],[0,0], color='black')
                if counter==1:
                    plt.ylabel("DF/F: \nBarrel Cortex (left) ", fontsize=10)
                counter+=1
            
            v_min = min(np.min(img[0]),np.min(img[1]),np.min(img[2]))
            v_max = max(np.max(img[0]),np.max(img[1]),np.max(img[2]))
            
            for k in range(3):
                ax=plt.subplot(3,3,3+1+k)
                plt.imshow(img[k], vmin=v_min,vmax=v_max)
                n_spikes = map_name[0][-9:-4]
                plt.title(plot_string+ "\n#spks: "+str(n_spikes))
                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)
                if k==0:
                    ax.get_yaxis().set_visible(True)
                    plt.ylabel("Relative DF/F", fontsize=10)

            
            plt.suptitle(file_name + "    Unit: " + str(unit), fontsize = 10)
            #plt.get_current_fig_manager().window.showMaximized()
            figManager = plt.get_current_fig_manager()
            figManager.window.showMaximized()
            plt.savefig(file_dir + file_name + '/unit_' +str(unit).zfill(2)+'.png')
            plt.close()
            #plt.show()
                
                
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
