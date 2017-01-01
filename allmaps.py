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
dropbox_dirs = []
file_names = []

main_dir = '/media/cat/12TB/in_vivo/tim/dongsheng/'
dropbox_dirs.append('/home/cat/Dropbox/tim/neuron_report/2015-5-6/')

#########**************** LOAD ALL FILES DATA IN EXPERIMENT DIRECTORY ****************

file_dirs.append('/media/cat/12TB/in_vivo/tim/dongsheng/2015-5-6/')  #*************************************
file_names.append([
'2015-5-6-1',       #cortex
'2015-5-6-2',       #cortex
'2015-5-6-3',      #subcortex
'2015-5-6-4',      #subcortex
'2015-5-6-5',      #subcortex
'2015-5-6-6'       #subcortex
])


plot_string='all'

window = 3

#************************************ Plot Max / Min Locations overlayed onto figures 
if False:
    dir_counter = 0
    file_dir=file_dirs[0]

    locations = ['barrel','allcortex']
    for location in locations: 
        max_pixel = []
        min_pixel = []

        for file_name in file_names[dir_counter]:
            #Load Quality control
            if False:
                temp_name = file_dir.replace('/media/cat/12TB/in_vivo/tim/','')
                dir_file_name = main_dir+"dongsheng_quality_control/"+temp_name+'/'+file_name+'/'
                files = os.listdir(dir_file_name)
                qc = np.zeros(100,dtype=np.int32)
                for f1 in files:
                    qc[int(f1[-6:-4])] = 1  #Save good units
                
            #Load depth
            depth = str(np.loadtxt(file_dir+file_name+'/depth.txt', dtype=str))
            print depth
                    
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
            
            for unit in units:
                traces=[]
                #if qc[int(unit)] == 0: continue 
                
                #Load max and min pixel locations:
                tc_name = glob.glob(file_dir + file_name+'/time_course_data_'+file_name+'_'+plot_string+'_'+str(window)+'sec_window_unit'+str(unit).zfill(2)+'*')
                if len(tc_name)==0: continue
                #print tc_name
                with open(tc_name[0], "r") as f:
                    data = csv.reader(f)
                    counter=0
                    for row in data:
                        traces.append(row)
                traces_temp = []
                indexes_temp = []
                min_pixel_temp=[]
                max_pixel_temp=[]
                for q in range(2,32,1):
                    traces_temp.append([float(z) for z in traces[q]])
                for q in range(32,62,2):
                    max_pixel_temp.append([int(z) for z in traces[q]])
                for q in range(33,62,2):
                    min_pixel_temp.append([int(z) for z in traces[q]])

                if location=='barrel':
                    min_pixel.append(min_pixel_temp[4])
                    max_pixel.append(max_pixel_temp[4])
                    
                elif location == 'allcortex':
                    if (np.isnan(traces_temp[28]).any()) or (np.isnan(traces_temp[29]).any()):         #Search for min
                        search_min = np.min(traces_temp[0])
                        temp_index = 0
                        for kk in range(0,30,2):
                            if np.min(traces_temp[kk])> search_min:
                                search_min=np.min(traces_temp[kk])
                                temp_index=kk/2
                        min_pixel.append(min_pixel_temp[temp_index])
                        #Search for max
                        search_max = np.max(traces_temp[0])
                        temp_index = 0
                        for kk in range(1,30,2):
                            if np.max(traces_temp[kk])> search_max:
                                search_max=np.min(traces_temp[kk])
                                temp_index=(kk-1)/2
                        max_pixel.append(min_pixel_temp[temp_index])
                    else:
                        min_pixel.append(min_pixel_temp[14])
                        max_pixel.append(max_pixel_temp[14])
                                    
        print "Number of max_pixels: ", len(max_pixel)
        map_name = glob.glob(file_dir + file_name+'/'+file_name+'_maxmap_unit_'+str(units[1]).zfill(3)+'*.npy')

        if len(map_name)==0:  #Check the _all condition which is the latest version of the maxmaps
            super_img = np.zeros((256,256),dtype=np.float32)
            for uu in range(len(units)):
                #if qc[uu] == 0: continue 

                map_name = glob.glob(file_dir + file_name+'/'+file_name+'_all_maxmap_unit_'+str(uu).zfill(2)+'*.npy')
                temp_img=np.load(map_name[0])

                v_min = np.min(temp_img)
                v_max = np.max(temp_img)
                temp_img = (temp_img - v_min)/(v_max-v_min)
                super_img+=temp_img

        img=super_img

        #Find max/min values of map to normalize data;
        img=np.load(map_name[0])
        v_min =np.min(img)
        v_max = 0
        for k in range(len(img)):
            #if (k>(len(img)*.35)) and (k<(len(img)*.45)): continue
            #temp = max(np.max(img[k,0:len(img)*.4]),np.max(img[k,len(img)*.6:]))
            #temp = np.max(img[k,0:len(img)*.4])
            temp = np.max(img[k])
            if v_max<temp: v_max = temp

        #********** PLOT Maxima ***********************
        pixel_width = 5

        ax =plt.subplot(1,2,1)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        img_max = np.float32(img.copy())*(-1)+1
        v_min = np.min(img_max)
        v_max = np.max(img_max)

        for u in range(len(max_pixel)):
            #for k in range(len(max_pixel[u])):
                for p in range(pixel_width):
                    for r in range(pixel_width):
                        img_max[max_pixel[u][0]-1+p,max_pixel[u][1]-1+r]= -100.

        v1a = masked_array(img_max,img_max<-50.)
        plt.imshow(v1a, vmin=v_min, vmax=v_max, alpha=1, cmap="Greys")
        v1b = masked_array(img_max,img_max>-50.)
        plt.imshow(v1b, vmin=-1000, vmax=-100, alpha=1, cmap="Reds")
        
        #********* PLOT Minima **********************
        ax =plt.subplot(1,2,2)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        img_min = np.float32(img.copy())*(-1)+1
        v_min = np.min(img_min)
        v_max = np.max(img_min)

        for u in range(len(min_pixel)):
            #for k in range(len(min_pixel[u])):
                for p in range(pixel_width):
                    for r in range(pixel_width):
                        img_min[min_pixel[u][0]-1+p,min_pixel[u][1]-1+r]= -100

        v1a = masked_array(img_min,img_min<-50.)
        plt.imshow(v1a, vmin=v_min, vmax=v_max, alpha=1, cmap="Greys")
        v1b = masked_array(img_min,img_min>-50.)
        plt.imshow(v1b, vmin=-1000, vmax=-100, alpha=1, cmap="Blues")
                        
        n_spikes = map_name[0][-9:-4]
        title_string = file_names[0][0]
        for i in range(1,len(file_names[0]),1):
            print file_names[0][i]
            title_string = title_string + "   " + file_names[0][i]
        print title_string
        title_string = textwrap.fill(title_string, 90)
        plt.suptitle(title_string + "\n "+depth + " max/min pixel locations", fontsize=10)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        #plt.ylabel("Normalized DF/F", fontsize=10)
        #plt.get_current_fig_manager().window.showMaximized()

        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        temp_name = file_dir.replace('/media/cat/12TB/in_vivo/tim/','')
        plt.savefig(file_dir + '/'+temp_name[:-1]+'_maxminpixel_'+depth+"_"+location+'.png')
        plt.close()
        #plt.show()
        
#*********************************** SAVE MAX MAPS AS FIGURES FOR QUALITY CONTROL *****************************
if True: 
    dir_counter = 0
    for file_dir, dropbox_dir in zip(file_dirs, dropbox_dirs):
    #for file_dir in file_dirs:
        print file_dir
        print file_names[dir_counter]
        #quit()

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
            
            units=np.array(units)
            indexes = np.argsort(units)
            units=units[indexes]
            channels=np.array(channels)[indexes]
            ptps=np.array(ptps)[indexes]
            
            #fig = plt.
            for unit in units:
                ax =plt.subplot(111)
                map_name = glob.glob(file_dir + file_name+'/'+file_name+'_maxmap_unit_'+str(unit).zfill(2)+'*.npy')

                if len(map_name)==0:  #Check the _all condition which is the latest version of the maxmaps
                    map_name = glob.glob(file_dir + file_name+'/'+file_name+'_all_maxmap_unit_'+str(unit).zfill(2)+'*.npy')
                    if len(map_name)==0: continue


                #print map_name[0]
                img=np.load(map_name[0])
                v_max = 0
                for k in range(len(img)):
                    #if (k>(len(img)*.35)) and (k<(len(img)*.45)): continue
                    #temp = max(np.max(img[k,0:len(img)*.4]),np.max(img[k,len(img)*.6:]))
                    #temp = np.max(img[k,0:len(img)*.4])
                    temp = np.max(img[k])
                    if v_max<temp: v_max = temp
        
                plt.imshow(img, vmin=np.min(img), vmax=v_max)
                n_spikes = map_name[0][-9:-4]
                plt.title(file_name+"  unit: " + str(unit)+"  ch: "+str(channels[unit])+"  #spks: "+str(n_spikes)+"\n v_min: "+str(round(np.min(img),4)*100)+"% v_max: "+str(round(v_max,4)*100)+"%", fontsize=15)
                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)
                #ax.get_yaxis().set_visible(True)
                plt.ylabel("Normalized DF/F", fontsize=10)
                #plt.get_current_fig_manager().window.showMaximized()
                figManager = plt.get_current_fig_manager()
                figManager.window.showMaximized()
                plt.savefig(dropbox_dir + file_name + '/unit_' +str(unit).zfill(2)+'.png')
                plt.close()
                #plt.show()


