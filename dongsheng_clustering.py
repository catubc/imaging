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
import sys
import os.path
import glob #Wildcard searching in dir listing
from pylab import *
import struct, array, csv
import pandas as pd
import matplotlib.mlab as mlab
import itertools
import cPickle as pkl
from scipy import interpolate
import dill  

from mpl_toolkits.mplot3d import Axes3D


main_dir = '/media/cat/12TB/in_vivo/tim/dongsheng/'

file_dirs = []
file_names = []

#########**************** LOAD ALL FILES DATA IN EXPERIMENT DIRECTORY ****************
file_dirs.append('/media/cat/12TB/in_vivo/tim/dongsheng/2015-5-6/')  #*************************************
file_names.append([
'2015-5-6-1',       #cortex
'2015-5-6-2',       #cortex
'2015-5-6-3',      #subcortex
'2015-5-6-4',      #subcortex
'2015-5-6-5',      #subcortex
'2015-5-6-6'       #subcortex

#'2015-5-6-7-FLstim',
#'2015-5-6-8-hlstim',
#'2015-5-6-9-v1',
#'2015-5-6-10-A1',
])

file_dirs.append('/media/cat/12TB/in_vivo/tim/dongsheng/2015-7-22/')  #*************************************
file_names.append([
'2015-7-22-1',     #cortex
'2015-7-22-2',     #cortex
'2015-7-22-3',     #cortex
'2015-7-22-4',     #subcortex
'2015-7-22-5',     #subcortex
'2015-7-22-6',     #subcortex
'2015-7-22-7',     #subcortex
'2015-7-22-8'      #subcortex

#'2015-7-22-9-W1', 
#'2015-7-22-10-w2', 
#'2015-7-22-11-w3',
#'2015-7-22-12-v1',
#'2015-7-22-13-V2',
#'2015-7-22-14-FL',
#'2015-7-22-15-FL',
#'2015-7-22-16-HL',
#'2015-7-22-17-HL',
#'2015-7-22-18-A1',
#'2015-7-22-19-A2'
])

file_dirs.append('/media/cat/12TB/in_vivo/tim/dongsheng/2015-7-23/')  #*************************************
file_names.append([
'2015-7-23-1',     #cortex
'2015-7-23-2',     #subcortex
'2015-7-23-3',     #subcortex  
'2015-7-23-16'     #subcortex

#'2015-7-23-5-w1',
#'2015-7-23-6-w2',
#'2015-7-23-7-v1',
#'2015-7-23-8-v2',
#'2015-7-23-9-a1',
#'2015-7-23-10-A2',
#'2015-7-23-11-FL',
#'2015-7-23-12-FL',
#'2015-7-23-13-HL',
#'2015-7-23-14-HL',
])

file_dirs.append('/media/cat/12TB/in_vivo/tim/dongsheng/2015-7-30/')  #*************************************
file_names.append([
'2015-7-30-4',
'2015-7-30-5', 
'2015-7-30-6', 
'2015-7-30-7'
])

file_dirs.append('/media/cat/12TB/in_vivo/tim/dongsheng/2015-11-18/') #*************************************
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


file_dirs.append('/media/cat/12TB/in_vivo/tim/dongsheng/2015-11-27/') #*************************************
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

file_dirs.append('/media/cat/12TB/in_vivo/tim/dongsheng/2015-12-1/') #*************************************
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

file_dirs.append('/media/cat/12TB/in_vivo/tim/dongsheng/2015-12-2/') #*************************************
file_names.append([
'2015-12-2-2-10electrodeincortex-iso1', 
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

file_dirs.append('/media/cat/12TB/in_vivo/tim/dongsheng/2015-12-11/') #*************************************
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

file_dirs.append('/media/cat/12TB/in_vivo/tim/dongsheng/2015-12-16/') #*************************************
file_names.append([
'2015-12-16-1-11electrodeincortex-iso1',
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


file_dirs.append('/media/cat/12TB/in_vivo/tim/dongsheng/2016-1-11/') #*************************************
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

file_dirs.append('/media/cat/12TB/in_vivo/tim/dongsheng/2016-1-14/') #*************************************
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

    

#****************** Cluster Single Area Time Traces
if True:
        
    n_areas = 1 #Number of cortical areas for which time courses computed

    window = 3  #Number of sec in computed time course window

    #area_names = [ 'hindlimb']
    area_names = ['barrel'] 
    maxmin= ['min','max']
    hemisphere=['right','left'] #Ok like this, index in time_course starts at 1 so need to start backwards for the mod function to work

    
    class STMTD(object):
        pass
    stmtd = STMTD()

    normalize = False

    stmtd.tc_data = []
    stmtd.area_maxormin = []
    stmtd.area_hemisphere = []

    stmtd.n_spikes = []
    stmtd.ptp = []
    stmtd.DF = []
    stmtd.depth_left = []
    stmtd.depth_right = []
    stmtd.state_left = []
    stmtd.state_right = []
    stmtd.mouse_index = []

    stmtd.channel_left = []
    stmtd.raster_left = []
    stmtd.template_left = []
    stmtd.area_index_left = []

    stmtd.file_index=[]
    stmtd.dir_index=[]
    exp_dirs=file_dirs
    sub_dirs=file_names
    
    load_traces(stmtd, n_areas, window, area_names, maxmin, hemisphere, file_dirs, file_names, main_dir)

    state='awake'
    depth='subcortical'

    n_spikes_thresh = 128 #min # of spikes
    ptp_thresh = 20 #min ptp amplitude in uV
    #minDF_thresh = 0.005 #min DF/F response required for inclusion ; 0.5 %
    minDF_thresh = 0.005 #min DF/F response required for inclusion ; 0.5 %
    
    cell_counter = 0
    good_traces=[]
    ave_trace = []
    marker = []
    xx = np.linspace(-3,3., 600)
    for k in range(len(stmtd.tc_data)):
        if (stmtd.n_spikes[k]>n_spikes_thresh) and (stmtd.ptp[k]>ptp_thresh) and (stmtd.DF[k]>minDF_thresh):
            #if (stmtd.depth_left[k]==depth) and (stmtd.state_left[k]==state):
            if True:
                max_val = max(np.max(stmtd.tc_data[k][200:400]),-np.min(stmtd.tc_data[k][200:600]))
                if stmtd.tc_data[k][0]/max_val>.5: os.remove(stmtd.dir_index[k]+'/'+stmtd.file_index[k])
                elif stmtd.tc_data[k][0]/max_val<-.5: os.remove(stmtd.dir_index[k]+'/'+stmtd.file_index[k])
                elif stmtd.tc_data[k][510]/max_val>.5: os.remove(stmtd.dir_index[k]+'/'+stmtd.file_index[k])
                elif stmtd.tc_data[k][565]/max_val>.5: os.remove(stmtd.dir_index[k]+'/'+stmtd.file_index[k])
                elif stmtd.tc_data[k][595]/max_val>.3: os.remove(stmtd.dir_index[k]+'/'+stmtd.file_index[k])
                
                else: 
                    good_traces.append(np.array(stmtd.tc_data[k])/max_val)
                    #good_traces.append(np.array(stmtd.tc_data[k])) #NON_Normalized traces
                    plt.plot(xx, np.array(stmtd.tc_data[k])/max_val, color='black', linewidth=1, alpha=0.1)
                    ave_trace.append(np.array(stmtd.tc_data[k])/max_val)
                    if stmtd.depth_left[k]=='subcortical': marker.append('o')
                    else: marker.append('^')
                    cell_counter+=1
                    
    plt.plot(xx, np.average(ave_trace,axis=0), color='black', linewidth=5, alpha=1.0)
    plt.plot([0,0],[-1,1], color='black')
    plt.plot([-3,3],[0,0], color='black')
    plt.ylim(-1,1)
    plt.xlim(-3,3)        
    
    print "# cells: ", cell_counter
    plt.title(depth+" "+state + " "+area_names[0]+ " #traces: "+str(len(ave_trace)))
    plt.show()
    
    if state=='anesthetized' and depth=='cortex': n_clusters = 1
    if state=='anesthetized' and depth=='subcortical': n_clusters = 4
    if state=='awake' and depth=='cortex': n_clusters = 2
    if state=='awake' and depth=='subcortical': n_clusters = 5
        
    cluster_traces(good_traces,stmtd, state, depth, n_clusters, marker)
    


quit()





#****************** MULTI-AREA ROI CLUSTERING ************
if True:
    
    dir_counter = 0
    for file_dir in file_dirs:
        print file_dir
        print file_names[dir_counter]
        for file_name in file_names[dir_counter]:

            #Load units from .csv file name; NB: This may lead to duplication as resaved .csv's for Dongsheng
            files = os.listdir(file_dir+file_name)
            temp_names = []
            for file_ in files:
                if ("unit_" in file_):
                    if ('map' in file_) or ('npy' in file_) or ('png' in file_):
                        pass
                    else:
                        temp_names.append(file_)

            #Save individual unit names, channels, ptps
            units = []
            channels = []
            ptps = []
            for file_ in temp_names:
                temp_unit_id = int(file_[5:7])

                #Skip units that failed qc
                if False:
                    if qc[temp_unit_id]==0: continue

                #Load good units
                units.append(temp_unit_id)
                channels.append(int(file_[16:18]))
                ptps.append(int(file_[23:26]))

            print len(units)
            quit()
            
            


