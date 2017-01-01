from tsf_ptcs_classes import *
from distributions import *
from sequential_firing import *

from scipy import stats
from scipy import signal
from numpy.ma import masked_array

import numpy as np
import time, math
import sys
import os.path
import multiprocessing as mp
import glob #Wildcard searching in dir listing
import textwrap

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
main_dir = '/media/cat/12TB/in_vivo/tim/'


#******************** LOAD DATA **********************
file_dirs.append('/media/cat/12TB/in_vivo/tim/2015-11-18/') #*************************************
file_names.append([
#'2015-11-18-1-9electrodein', 
#'2015-11-18-2-9electrodein', 
#'2015-11-18-3-9electrodein',
#'2015-11-18-4-9electrodein-iso0.5',
#'2015-11-18-5-9electrodein-iso0.5',
#'2015-11-18-6-9electrodein-iso0.5',
#'2015-11-18-7-9electrodein-iso0',
#'2015-11-18-8-9electrodein-iso0',
#'2015-11-18-9-9electrodein-iso0',
#'2015-11-18-10-allelectrodein-iso0',
#'2015-11-18-11-allelectrodein-iso0',
#'2015-11-18-12-allelectrodein-iso0',
#'2015-11-18-13-deep-iso0',
'2015-11-18-14-deep-iso0', 
#'2015-11-18-15-deep-iso0', 
])

plot_string='all'
window = 3

#************************************ Plot Max / Min Locations overlayed onto figures 
dir_counter = 0
file_dir=file_dirs[0]

n_pixels=256
#Load General mask (removes background)
generic_mask_file = []
generic_mask_file = file_dir + 'genericmask.txt'
if (os.path.exists(generic_mask_file)==True):
    generic_coords = np.loadtxt(generic_mask_file)

#Load Artifact mask 
artifact_coords = []
artifact_mask_file = file_dir + 'artifactmask.txt'
if (os.path.exists(artifact_mask_file)==True):
    artifact_coords = np.loadtxt(artifact_mask_file)

#Load Bregma coords for computing contralateral areas
bregma_mask_file = file_dir + 'bregmamask.txt'
bregma_coords = np.loadtxt(bregma_mask_file)
if len(bregma_coords)==2:
    bregma_coords = [bregma_coords]
bregma_coords_temp = []
for j in range(len(bregma_coords)):  #Remove centreline from all images
    for k in range(n_pixels):
        for l in range(7):
            bregma_coords_temp.append([k,min(n_pixels-1,bregma_coords[j][1]-3+l)])

generic_mask_indexes=np.zeros((n_pixels,n_pixels))
for i in range(len(generic_coords)):
    generic_mask_indexes[generic_coords[i][0]][generic_coords[i][1]] = True
for i in range(len(bregma_coords_temp)):
    generic_mask_indexes[bregma_coords_temp[i][0]][bregma_coords_temp[i][1]] = True       
for i in range(len(artifact_coords)):
    generic_mask_indexes[artifact_coords[i][0]][artifact_coords[i][1]] = True     


for file_name in file_names[dir_counter]:
    #Load Quality control
    temp_name = file_dir.replace('/media/cat/12TB/in_vivo/tim/','')
    dir_file_name = main_dir+"dongsheng_quality_control/"+temp_name+'/'+file_name+'/'
    files = os.listdir(dir_file_name)
    qc = np.zeros(100,dtype=np.int32)
    for f1 in files:
        qc[int(f1[-6:-4])] = 1  #Save good units

    #Load depth
    #depth = str(np.loadtxt(file_dir+file_name+'/depth.txt', dtype=str))
    #print depth
            
    #Load unit names
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
        #print file_
        units.append(int(file_[5:7]))
        channels.append(int(file_[16:18]))
        ptps.append(int(file_[23:26]))
    
    units=np.array(units)
    indexes = np.argsort(units)
    units=units[indexes]
    channels=np.array(channels)[indexes]
    ptps=np.array(ptps)[indexes]

    max_pixel = []
    min_pixel = []
    
    for unit in units:
        traces=[]
        if qc[int(unit)] == 0: continue 
        if unit!=2: continue
        #Load max and min pixel locations:
        #tc_name = glob.glob(file_dir + file_name+'/time_course_data_'+file_name+'_unit'+str(unit).zfill(2)+'_ch'+str(channels[unit]))
        #print file_dir + file_name+'/time_course_data_'+file_name+'_'+str(unit).zfill(2)+'*'
        tc_name = glob.glob(file_dir + file_name+'/time_course_data_'+file_name+'_'+plot_string+'_'+str(window)+'sec_window_unit'+str(unit).zfill(2)+'*')
        print tc_name
        
        if len(tc_name)==0: continue
        #print tc_name
        with open(tc_name[0], "r") as f:
            data = csv.reader(f)
            counter=0
            for row in data:
                traces.append(row)

        traces_temp=[]
        indexes_temp = []
        min_pixel_temp=[]
        max_pixel_temp=[]
        
        for q in range(2,32,1):
            traces_temp.append([float(z) for z in traces[q]])
        for q in range(32,62,2):
            print q
            max_pixel_temp.append([int(z) for z in traces[q]])
        for q in range(33,62,2):
            min_pixel_temp.append([int(z) for z in traces[q]])

        location='barrel'
        if location=='barrel':
            min_pixel=min_pixel_temp[4]
            max_pixel=max_pixel_temp[4]
        print "max_pixel: ", max_pixel
        
        traces_temp=np.array(traces_temp)
        if np.max(abs(traces_temp[8:10]))<0.0075: continue
        
        if np.max(abs(traces_temp[8]))>np.max(abs(traces_temp[9])): 
            traces_plot=traces_temp[8]
            norm_val='min'
        else:
            traces_plot=traces_temp[9]
            norm_val='max'
        
        
        color='red'
        #plt.plot(traces_temp[6], color='cyan')
        #plt.plot(traces_temp[7], color='pink')
        #plt.plot(traces_temp[8], color='blue')
        xx = np.linspace(-3.0,3.0,len(traces_temp[0]))
        ax=plt.subplot(1,2,1)
        plt.plot(xx, traces_plot, color=color, linewidth=5)
        plt.ylim(-0.02,0.06)
        plt.plot([-3,3],[0,0],color='black',linewidth=3)
        plt.plot([0,0],[-1,1],color='black',linewidth=3)

        ax=plt.subplot(1,2,2)
        plt.plot(xx, traces_plot, color=color, linewidth=5)
        ax.axis('off')
        plt.ylim(-0.02,0.06)
        plt.plot([-3,3],[0,0],color='black',linewidth=3)
        plt.plot([0,0],[-1,1],color='black',linewidth=3)

        plt.suptitle(file_name+ "   unit: "+str(unit))
        temp_name = file_dir.replace('/media/cat/12TB/in_vivo/tim/','')
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        plt.savefig(file_dir + '/'+file_name+'_unit_'+ str(unit)+'_framestraces.png',dpi=300)

        #plt.show()
        #for q in range(32,62,2):
            #max_pixel_temp.append([int(z) for z in traces[q]])
        #for q in range(33,62,2):
            #min_pixel_temp.append([int(z) for z in traces[q]])

        #if location=='barrel':
            #min_pixel.append(min_pixel_temp[4])
            #max_pixel.append(max_pixel_temp[4])
                

        #map_name = glob.glob(file_dir + file_name+'/'+file_name+'_all_maxmap_unit_'+str(uu).zfill(2)+'*.npy')
        images_name = glob.glob(file_dir+file_name+'/img_avg_'+file_name+'_unit'+str(unit).zfill(2)+'*.npy')
        print file_dir+file_name+'/img_avg_'+file_name+'_unit'+str(unit).zfill(2)+'*.npy'
        print images_name
        
        images = np.load(images_name[0])
        #plt.close()
        #temp_array = np.ma.array(images[k], mask=generic_mask_indexes, fill_value =0, hard_mask = True)

        #plt.imshow(images[0])
        #plt.show()
        v_max = 0
        v_min = 0
        chunk_start=[0,60,120]
        chunk_end=[60,120,178]
        pixel_width=20
        
        if norm_val=='max':
            blank_pixel=max_pixel
        else:
            blank_pixel=min_pixel
                    
        for k in range(0,len(images),1):
            temp_array = np.ma.array(images[k], mask=generic_mask_indexes, fill_value =0, hard_mask = True)
            if np.nanmax(temp_array)>v_max: v_max =np.nanmax(temp_array)
            if np.nanmin(temp_array)<v_min: v_min =np.nanmin(temp_array)
        v_max = v_max
        #v_max = 0.01
        v_min = v_min
            
        for chunk1, chunk2 in zip(chunk_start,chunk_end):
            for k in range(chunk1,chunk2,1):
                ax=plt.subplot(6,10,k-chunk1+1)
                temp_array = np.ma.array(images[k], mask=generic_mask_indexes, fill_value =0, hard_mask = True)
                
                for p in range(pixel_width):
                    for r in range(pixel_width):
                        temp_array[blank_pixel[0]-pixel_width/2+p,blank_pixel[1]-pixel_width/2+r]= v_max #(v_max+v_min)/2.
                plt.imshow(temp_array, vmin=v_min, vmax=v_max)
                plt.title(str(round(float(k)/30.-3,2))+" sec.", fontsize=5)

                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)
            plt.suptitle(file_name+ "  unit: "+str(unit)+ "   DF/F: vmin: "+str(round(v_min*1E2,1))+ "  vmax: "+str(round(v_max*1E2,1)))

            #temp_name = file_dir.replace('/media/cat/12TB/in_vivo/tim/','')
            figManager = plt.get_current_fig_manager()
            figManager.window.showMaximized()
            plt.savefig(file_dir + '/'+file_name+'_unit_'+str(unit)+ 'frames_'+str(chunk1/60)+'.png',dpi=300)
            #plt.show()

        #quit()

