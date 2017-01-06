import sys
sys.path.insert(0, '/home/cat/code/ephys')
from tsf_ptcs_classes import *
from distributions import *
#from sequential_firing import *

from scipy import stats
from scipy import signal

import numpy as np
import time, math
import os.path
import multiprocessing as mp
import re

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


main_dir = '/media/cat/8TB/in_vivo/tim/dongsheng/'

file_dirs = []
file_names = []

#############**************** LOAD ALL FILES DATA IN EXPERIMENT DIRECTORY ****************
#file_dirs.append('/media/cat/8TB/in_vivo/tim/dongsheng/2015-5-6/')  #*************************************
#file_names.append([
#'2015-5-6-1',       #cortex
#'2015-5-6-2',       #cortex
#'2015-5-6-3',      #subcortex      #DO NOT USE THIS ANIMAL - SMALLER CRANIUM CAUSES WEIRD OFFSCREEN ARTIFACTS
#'2015-5-6-4',      #subcortex
#'2015-5-6-5',      #subcortex
#'2015-5-6-6'       #subcortex

##'2015-5-6-7-FLstim',
##'2015-5-6-8-hlstim',
##'2015-5-6-9-v1',
##'2015-5-6-10-A1',
#])

#file_dirs.append('/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/')  #*************************************
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

#file_dirs.append('/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-23/')  #*************************************
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

#file_dirs.append('/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-30/')  #*************************************
#file_names.append([
#'2015-7-30-4', #subcortex
#'2015-7-30-5', #subcortex
#'2015-7-30-6', #subcortex
#'2015-7-30-7'  #subcortex
#])

#file_dirs.append('/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/') #*************************************
#file_names.append([
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
#'2015-11-18-14-deep-iso0', 
#'2015-11-18-15-deep-iso0', 
#'2015-11-18-16-deep-iso0', 
#'2015-11-18-17-deep-iso0', 
#'2015-11-18-18-deep-iso0'
#])


file_dirs.append('/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/') #*************************************
file_names.append([
#'2015-11-27-1-10electrodein-iso1.0',
#'2015-11-27-2-10electrodein-iso1', 
#'2015-11-27-4-10electrodein-iso0', 
'2015-11-27-5-10electrodein-iso0',
#'2015-11-27-7-16electrodein-iso1',
#'2015-11-27-8-16electrodein-iso1',
#'2015-11-27-10-16electrodein-iso0',
#'2015-11-27-11-16electrodein-iso0',
#'2015-11-27-13-deep-iso1',
#'2015-11-27-14-deep-iso1',
#'2015-11-27-16-deep-iso0',

###'2015-11-27-2-10electrodein-iso1_mua', 
###'2015-11-27-5-10electrodein-iso0_mua',
###'2015-11-27-10-16electrodein-iso0_mua',
###'2015-11-27-14-deep-iso1_mua',
###'2015-11-27-16-deep-iso0_mua'

###'2015-11-27-2-10electrodein-iso1_lfp', 
###'2015-11-27-5-10electrodein-iso0_lfp',
###'2015-11-27-14-deep-iso1_lfp',
###'2015-11-27-16-deep-iso0_lfp'
])

#file_dirs.append('/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/') #*************************************
#file_names.append([
#'2015-12-1-1-10electrodeiniso1',
#'2015-12-1-2-10electrodeiniso1',
#'2015-12-1-3-10electrodeiniso1',
#'2015-12-1-5-10electrodeiniso0',
#'2015-12-1-6-10electrodeiniso0',
#'2015-12-1-8-allelectrodeiniso0.8', 
#'2015-12-1-9-allelectrodeiniso0.8', 
#'2015-12-1-11-allelectrodeiniso0', 
#'2015-12-1-12-allelectrodeiniso0',
#'2015-12-1-14-5electrodeinthalamus-iso0.8',
#'2015-12-1-15-5electrodeinthalamus-iso0.8',
#'2015-12-1-18-5electrodeinthalamus-iso0',   #SINGLE
#'2015-12-1-19-allelectrodeinthalamus-iso0.8',  #SINGLE
#'2015-12-1-20-allelectrodeinthalamus-iso0.8',
#'2015-12-1-22-allelectrodeinthalamus-iso0',
#'2015-12-1-23-allelectrodeinthalamus-iso0', 
#'2015-12-1-24-allelectrodeinthalamus-iso0'

###'2015-12-1-3-10electrodeiniso1_mua',
###'2015-12-1-6-10electrodeiniso0_mua',
###'2015-12-1-20-allelectrodeinthalamus-iso0.8_mua',
###'2015-12-1-23-allelectrodeinthalamus-iso0_mua', 

###'2015-12-1-3-10electrodeiniso1_lfp',
###'2015-12-1-6-10electrodeiniso0_lfp',
###'2015-12-1-20-allelectrodeinthalamus-iso0.8_lfp',
###'2015-12-1-23-allelectrodeinthalamus-iso0_lfp', 
#])

#file_dirs.append('/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/') #*************************************
#file_names.append([
#'2015-12-2-2-10electrodeincortex-iso1', 
#'2015-12-2-3-10electrodeincortex-iso1',
#'2015-12-2-4-10electrodeincortex-iso1',
#'2015-12-2-6-10electrodeincortex-iso0',
#'2015-12-2-7-10electrodeincortex-iso0',
#'2015-12-2-9-allelectrodeincortex-iso1',  #SINGLE
#'2015-12-2-11-allelectrodeinthalamus-iso1',
#'2015-12-2-12-allelectrodeinthalamus-iso1',
#'2015-12-2-14-allelectrodeinthalamus-is0',
#'2015-12-2-15-allelectrodeinthalamus-is0',

##'2015-12-2-3-10electrodeincortex-iso1_mua',
##'2015-12-2-6-10electrodeincortex-iso0_mua',
##'2015-12-2-12-allelectrodeinthalamus-iso1_mua',
#'2015-12-2-14-allelectrodeinthalamus-is0_mua',

##'2015-12-2-3-10electrodeincortex-iso1_lfp',
##'2015-12-2-6-10electrodeincortex-iso0_lfp',
##'2015-12-2-12-allelectrodeinthalamus-iso1_lfp',
#'2015-12-2-14-allelectrodeinthalamus-is0_lfp'
#])

#file_dirs.append('/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-11/') #*************************************
#file_names.append([
#'2015-12-11-1-10electincortex-iso1.0',
#'2015-12-11-2-10electincortex-iso1',
#'2015-12-11-3-10electincortex-iso1',
#'2015-12-11-5-10electincortex-iso0',
#'2015-12-11-7-10electincortex-iso1',
#'2015-12-11-8-10electincortex-iso1',
#'2015-12-11-9-allelectincortex',
#'2015-12-11-10-5electinthalamusorcpu',
#'2015-12-11-11-allelectinthalamus-iso1',
#'2015-12-11-12-allelectinthalamus-iso1',
#'2015-12-11-13-allelectinthalamus-iso1',
#'2015-12-11-15-allelectinthalamus-iso0',
#'2015-12-11-16-allelectinthalamus-iso0',
#'2015-12-11-17-allelectinthalamus-iso0',
#'2015-12-11-18-allelectinthalamus-iso1'

##'2015-12-11-3-10electincortex-iso1_mua',
#'2015-12-11-5-10electincortex-iso0_mua',
##'2015-12-11-12-allelectinthalamus-iso1_mua',
##'2015-12-11-16-allelectinthalamus-iso0_mua',

#'2015-12-11-3-10electincortex-iso1_lfp',
#'2015-12-11-5-10electincortex-iso0_lfp',
##'2015-12-11-12-allelectinthalamus-iso1_lfp',
#'2015-12-11-16-allelectinthalamus-iso0_lfp',
#])

#file_dirs.append('/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-16/') #*************************************
#file_names.append([
#'2015-12-16-1-11electrodeincortex-iso1',
#'2015-12-16-2-11electrodeincortex-iso1',
#'2015-12-16-4-11electrodeincortex-iso0',
#'2015-12-16-5-11electrodeincortex-iso0',
#'2015-12-16-7-11electrodeincortex-iso1',
#'2015-12-16-8-4electrodeincpu-iso1',
#'2015-12-16-9-allelectrodeinZI-iso1',
#'2015-12-16-10-allelectrodeinZI-iso1',
#'2015-12-16-11-allelectrodeinZI-iso1',
#'2015-12-16-13-allelectrodeinZI-iso0',
#'2015-12-16-14-allelectrodeinZI-iso0',
#'2015-12-16-16-allelectrodeinZI-iso1',
#'2015-12-16-17-allelectrodeinZI-iso1'
#])

#file_dirs.append('/media/cat/8TB/in_vivo/tim/dongsheng/2016-1-11/') #*************************************
#file_names.append([
#'2016-1-11-1-10electrodeincortex-iso1.2',
#'2016-1-11-2-10electrodeincortex-iso1',
#'2016-1-11-3-10electrodeincortex-iso1',
#'2016-1-11-5-10electrodeincortex-iso0',
#'2016-1-11-6-10electrodeincortex-iso0',
#'2016-1-11-8-10electrodeincortex-iso1',
#'2016-1-11-9-allelectrodeincortex-iso1',
#'2016-1-11-10-13electrodeinthalamus-iso1',
#'2016-1-11-11-13electrodeinthalamus-iso1',
#'2016-1-11-12-13electrodeinthalamus-iso1',
#'2016-1-11-14-13electrodeinthalamus-iso0',
#'2016-1-11-15-13electrodeinthalamus-iso0',
#'2016-1-11-17-13electrodeinthalamus-iso1',
#'2016-1-11-18-13electrodeinthalamus-iso1'
#])

#file_dirs.append('/media/cat/8TB/in_vivo/tim/dongsheng/2016-1-14/') #*************************************
#file_names.append([
#'2016-1-14-1-10electrodein cortex-iso1',
#'2016-1-14-2-10electrodein cortex-iso1',
#'2016-1-14-3-10electrodein cortex-iso1',
#'2016-1-14-5-10electrodein cortex-iso0',
#'2016-1-14-6-10electrodein cortex-iso0',
#'2016-1-14-7-10electrodein cortex-iso0',
#'2016-1-14-9-10electrodein cortex-iso1',
#'2016-1-14-10-allelectrodein cortex-iso1',
#'2016-1-14-11-allelectrodeinthalamus-iso1',
#'2016-1-14-12-allelectrodeinthalamus-iso1',
#'2016-1-14-13-allelectrodeinthalamus-iso1',
#'2016-1-14-15-allelectrodeinthalamus-iso0',
#'2016-1-14-16-allelectrodeinthalamus-iso0',
#'2016-1-14-17-allelectrodeinthalamus-iso0'
#])

window = 3 #window to search in seconds
n_procs=14

area_names = ['hindlimb', 'forelimb', 'barrel','retrosplenial','visual', 'motor', 'pta', 'acc'] 
#area_names = ['hindlimb']

sides = ['left','right']

stm_types = ["all"]  # ["all", "burst", "1sec"]    #Do "all" last to ensure that the time courses are saved also

#area_names = [ 'area1']
#sides = ['left']

#STA maps - False units
sta_maps = True
overwrite = True                   #Compute stimulus maps 
overwrite_time_courses = False      #Flag to bipass existing img_avg generation and just compute time courses
overwrite_shifted=False             #Flag for overwriting rotated/shifted .npy image files

#Stimulus maps: define cortical areas using stimulus recordings
stimulus_maps = False

#Static maps - SUA and LFP
sua_static_maps = False             #make maxmap and minmap data files (required Clustering shapes of max value pixel)
lfp_static_maps = False

#Specificity maps - SUA, MUA and LFP
specificity_maps_over_depth = False
specificity_maps_over_time = False
specificity_maps_each_channel = False
specificity_maps_overlap_adjacent_channel = False
specificity_maps_lfp = False

#Make amax maps (matthieu method) for demarkating cortical areas
max_maps = False

#LFP Correlation maps - NOT USED ANYMORE; see lfp_static_maps which 
#lfp_correlation_maps = False 

max_map_list=[]
min_map_list=[]
max_map_nspikes=[]
min_map_nspikes=[]
max_map_channel=[]
min_map_channel=[]
max_map_ptp=[]
min_map_ptp=[]
max_map_depth=[]
min_map_depth=[]
max_map_state=[]
min_map_state=[]
max_map_unit=[]

dir_counter = 0
for file_dir in file_dirs:
    print file_dir
    print file_names[dir_counter]
    for file_name in file_names[dir_counter]:

        ##Load quality control file; exclude units that did not make QC!
        #if False:
            #temp_fname = file_dir.replace(main_dir, '')
            #dir_file_name = main_dir+"dongsheng_quality_control/"+temp_fname+file_name+'/'
            #files = os.listdir(dir_file_name)
            #qc = np.zeros(100,dtype=np.int32)
            #for f1 in files:
                ##print f1, f1[-6:-4]
                #qc[int(f1[-6:-4])] = 1  #S

        #Load units from .csv file name; NB: This may lead to duplication as resaved .csv's for Dongsheng
        files = os.listdir(file_dir+file_name)
        temp_names = []
        for file_ in files:
            if ("unit_" in file_):
                if ('map' in file_) or ('npy' in file_) or ('png' in file_) or ('parcel' in file_):
                    pass
                else:
                    temp_names.append(file_)

        #Save individual unit names, channels, ptps
        units = []
        channels = []
        ptps = []
        for file_ in temp_names:
            units.append(int(file_[5:7]))
            channels.append(int(file_[16:18]))
            ptps.append(int(file_[23:26]))

        #Load depth of recording and anesthetic state:
        infile = file_dir+file_name+"/state.txt"
        with open(infile, "r") as f:
            for line in f:
                state = line.rstrip()
        infile = file_dir+file_name+"/depth.txt"
        with open(infile, "r") as f:
            for line in f:
                depth = line.rstrip()


        #*********** Compute max maps from all cell pre-computed static maps
        if max_maps or specificity_maps_over_depth or specificity_maps_lfp or specificity_maps_over_time \
        or specificity_maps_each_channel or specificity_maps_overlap_adjacent_channel:
            temp = Load_max_maps(ptps, file_dir, file_name) #Returns: max_maps, min_maps, max_maps_nspikes, min_maps_nspikes
            print "Loading static maps: ", file_name
            max_map_list.append(temp[0])
            min_map_list.append(temp[1])
            max_map_nspikes.append(temp[2])
            min_map_nspikes.append(temp[3])
            max_map_ptp.append(temp[4])
            min_map_ptp.append(temp[5])
            max_map_depth.append(temp[6])
            min_map_depth.append(temp[7])
            max_map_state.append(temp[8])
            min_map_state.append(temp[9])
            max_map_channel.append(temp[10])
            min_map_channel.append(temp[11])
            max_map_unit.append(temp[12])
            
            continue

        #*********** Define Stimulus Triggered cortical areas
        if False:
            #Compute stimulus maps
            if False:
                Compute_stimulus_maps(main_dir, file_dir, area_names)
            
            if True:
            #Show all stimulus maps simultaneously; save maps; save master bregma/lambda locations
                roi_img = Show_master_stimulus_maps(main_dir, area_names)
                plt.imshow(roi_img)
                plt.plot([128,128],[0,256], linewidth=2, color='white')
                plt.show()
                
            quit()

        #*********** Define Cortical Areas
        if False:  
            Define_cortical_areas(file_dir, file_name, area_names, sides)
            continue
            
        #*********** Load raw images and rotate + align them
        if True:
            n_pixels = 256

            file_name_aligned = file_dir + file_name+'/'+file_name+'_images_aligned.npy'
            if (os.path.exists(file_name_aligned)==False):
                images_raw = Load_images(file_dir, file_name)
                images_rotated = Rotate_images(images_raw, file_dir, file_name, overwrite_shifted)
                images_aligned = Define_bregma_lambda(images_rotated, main_dir, file_dir, file_name)
            else:
                images_aligned = np.load(file_name_aligned)


            #*********** Load start and end times of imaging system and compute required parameters
            img_start, img_end, len_frame, img_rate, n_pixels, n_frames, img_times = Load_images_start_end(file_dir, file_name, images_aligned)
            print "Imaging rate: ", img_rate
            
            #Define_generic_mask(images_aligned, file_dir, file_name, window, img_rate, n_pixels)

        ##*********** Ablate Pixel Areas (i.e remove small artifacts)
        #if False:  Define_artifact(images_aligned, file_dir, file_name, window, n_pixels, img_rate)
        
        #************ Load generic mask
        if True:
            #Load General mask (removes background)
            generic_mask_file = main_dir + 'genericmask.txt'
            if (os.path.exists(generic_mask_file)==True):
                generic_coords = np.int16( np.loadtxt(generic_mask_file))
           
            generic_mask_indexes=np.int16(np.zeros((256,256)))
            for i in range(len(generic_coords)):
                generic_mask_indexes[generic_coords[i][0]][generic_coords[i][1]] = True
        #************ Load area borders for plotting
        if False:
            area_borders = []
            for area in area_names:
                for side in sides:
                    area_borders.append(np.load(main_dir + area+'_' + side + '_contour.npy'))
            
            #print len(area_borders)
            #quit()
        
        #*********** Compute Single Unit maps (vs. LFP maps)
        if sta_maps or sua_static_maps:
            print "No. units: ", len(units)
            print "OVERWRITE: ", overwrite

            #**************** LOOP OVER ALL UNITS IN EXPERIMENT **************************
            for i in range(len(units)):
            #for i in [0]:
                print "*Processing ", file_name, " unit: ", i+1, " of ", len(units), " unit name: ", units[i]
                unit=units[i]
                channel=channels[i]
                ptp=ptps[i]
                
                #PICK A SINGLE UNIT TO LOOK AT
                if unit not in [12]: continue
                
                #*********** Load sorted spikes for unit
                spikes = np.loadtxt(file_dir+file_name+'/unit_'+str(unit).zfill(2)+ '_channel_' + str(channel).zfill(2) + '_ptp_'+str(ptp).zfill(3)+'.csv')

                #Determine if any spikes in recording window to continue:
                temp0 = np.where(np.logical_and(spikes>=img_times[0]+window, spikes<=img_times[-1]-window))[0]
                spikes_in_window = spikes[temp0]
                if len(spikes_in_window) == 0: 
                    print "Zero spikes in window - skip unit"
                    continue
                n_spikes = len(spikes_in_window)  #NB: Save n_spikes value just for current epoch for continuous recordings

                #*********** Compute static maps +/- 1sec from spike rasters
                if sua_static_maps:  
                    Compute_static_maps_max(img_rate, window, n_procs, main_dir, file_dir, file_name, n_pixels, unit, channel, n_spikes, ptp)
                    continue   #If computing static maps: skip to next unit - do not reprocess data

                #*********** Compute spike triggered averages and remove baseline
                images_processed, spikes, plot_string = Compute_spike_triggered_average(unit, channel, spikes, window, img_rate, img_times, n_pixels, 
                                                            images_aligned, file_dir, file_name, n_procs, overwrite, stm_types)

                #*********** Compute spike-triggered maps on luminance-constrast feature map
                #if False:
                    #images_nobaseline = Remove_global_baseline(images_aligned)
                    ##images_contrast = Compute_luminance_contrast(images_nobaseline)
                    #images_areas_nobaseline, coords_save, generic_mask_indexes = Load_areas_whole_cortex(unit, channel, n_pixels, file_dir, file_name, images_nobaseline)
                    #images_decorrelated, tau = Compute_decorrelated_images(images_areas_nobaseline, file_dir, file_name, unit, channel, plot_string, window, spikes_in_window, n_procs, img_times, generic_mask_indexes)
                    #Max_plot, Min_plot, Max_pixel_value, Min_pixel_value, Max_index, Min_index = Search_max_min_single_area(images_decorrelated, img_rate, window, coords_save, n_procs, tau)
                    #Animate_images_decorrelated(unit, channel, window, images_decorrelated, file_dir, file_name, n_pixels, spikes, plot_string, n_procs, generic_mask_indexes,Max_pixel_value, Min_pixel_value)

                #*********** Load Areas, Search Max/Min, Save Time Courses
                if False:
                    #images_areas = Load_areas_and_mask(depth, unit, channel, n_pixels, main_dir, file_dir, file_name, images_processed, area_names, sides)

                    #*********** Search for Min and Max pixels
                    Max_plot, Min_plot, Max_pixel_value, Min_pixel_value, Max_index, Min_index, area_names, images_areas = \
                        Search_max_min(images_processed, file_dir, file_name, img_rate, window, n_procs, area_names, depth, sides, generic_mask_indexes)

                    #*********** Save time course .npy file and figures
                    print "Saving time course for unit : ", unit, " ..."
                    Save_time_course(unit, channel, spikes, Max_plot, Min_plot, Max_index, Min_index, window, len_frame, file_dir, file_name, area_names, sides, plot_string)

                    if False: 
                        Animate_images(unit, channel, window, img_rate, images_areas, file_dir, file_name, n_pixels, spikes, plot_string, n_procs, generic_mask_indexes, 
                        Max_plot, Min_plot, Max_pixel_value, Min_pixel_value, Max_index, Min_index, area_names, sides, depth)
                
                #*********** Generate ROI matrix and time courses
                if False:
                    
                    images_areas = Load_areas_and_mask(depth, unit, channel, n_pixels, main_dir, file_dir, file_name, images_processed, area_names, sides)
                    
                    average_areas = Average_roi(images_areas, img_rate, window, n_procs, area_names, sides)

                    Plot_matrix_maps(average_areas, file_dir, file_name, area_names, img_rate, unit, spikes, channel, ptp)


                print "...done."
                print " "
        #if lfp_correlation_maps:
            #n_pixels = 256

            #file_name_aligned = file_dir + file_name+'/'+file_name+'_images_aligned.npy'
            #if (os.path.exists(file_name_aligned)==False):
                #images_raw = Load_images(file_dir, file_name)
                #images_rotated = Rotate_images(images_raw, file_dir, file_name, overwrite_shifted)
                #images_aligned = Define_bregma_lambda(images_rotated, main_dir, file_dir, file_name)
            #else:
                #images_aligned = np.load(file_name_aligned)

            #lfp_img_correlation(file_dir, file_name, images_aligned,img_start,img_end,len_frame,img_rate,n_pixels,n_frames, img_times)
        
        if lfp_static_maps:
            n_pixels = 256

            file_name_aligned = file_dir + file_name+'/'+file_name+'_images_aligned.npy'
            if (os.path.exists(file_name_aligned)==False):
                images_raw = Load_images(file_dir, file_name)
                images_rotated = Rotate_images(images_raw, file_dir, file_name, overwrite_shifted)
                images_aligned = Define_bregma_lambda(images_rotated, main_dir, file_dir, file_name)
            else:
                images_aligned = np.load(file_name_aligned)
            
            Compute_lpf_static_maps(file_dir, file_name, images_aligned, img_start, img_end, len_frame, img_rate, n_pixels, n_frames, img_times)

    dir_counter+=1

if max_maps or specificity_maps_over_depth or specificity_maps_lfp or specificity_maps_over_time \
    or specificity_maps_each_channel or specificity_maps_overlap_adjacent_channel:
    pass
else:
    quit()


if specificity_maps_over_depth:

    #************* LOAD MASKS *****************
    #n_pixels = len(max_map_list[0][0])
    n_pixels = 256

    #Load General mask (removes background)
    generic_mask_file = []
    generic_mask_file = file_dir + 'genericmask.txt'
    if (os.path.exists(generic_mask_file)==True):
        generic_coords = np.int16(np.loadtxt(generic_mask_file))

    #Load Artifact mask 
    artifact_coords = []
    artifact_mask_file = file_dir + 'artifactmask.txt'
    if (os.path.exists(artifact_mask_file)==True):
        artifact_coords = np.int16(np.loadtxt(artifact_mask_file))

    #Load Bregma coords for computing contralateral areas
    bregma_mask_file = file_dir + 'bregmamask.txt'
    bregma_coords = np.int16(np.loadtxt(bregma_mask_file))
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
    #for i in range(len(bregma_coords_temp)):
    #    generic_mask_indexes[bregma_coords_temp[i][0]][bregma_coords_temp[i][1]] = True       
    for i in range(len(artifact_coords)):
        generic_mask_indexes[artifact_coords[i][0]][artifact_coords[i][1]] = True     
    
    midline_mask_n_pixels = 10
    for i in range(midline_mask_n_pixels):
        generic_mask_indexes[:,n_pixels/2+int(int(midline_mask_n_pixels)/2)-i]=True
    cortex_pixels = n_pixels**2-np.sum(generic_mask_indexes)

    min_ptp = 1
    if 'mua' in file_name:  #MUA should plot all mixed units regardless of amplitude
        min_ptp = 20
    min_spikes = 1
    
    ratios=[]
    
    depths = ['cortex','subcortical']
    states = ['anesthetized','awake']
    #map_list contains of list of all cells in each experiment: [recording, cell#]
    counter=1
    cell_coverage=[]

    #cell_list = [0,2,4,7,11,12,15,16,18,19,21,22,24,25,27,28,33,38,44,49,54,58]

    for depth in depths:
        for state in states:
            print "Processing: ", depth, state               
            ax = plt.subplot(2,2,counter)
            save_indexes=[]
            colors= []
            cell_counter=0
            #cell_coverage.append([])
            for channel_pick in np.arange(0,16,1):  #Loop over channels so can limit one map per channel
                channel_maps = []
                for i in range(len(max_map_list)):
                    #cell_break=False
                    for j in range(len(max_map_list[i])):       #Compute max map from all cells in each experiment
                        #print "unit, ch, i, j: ", max_map_unit[i][j], channel_pick, i, j
                        if (max_map_depth[i][j]==depth) and (max_map_state[i][j]==state) and (max_map_ptp[i][j]>min_ptp) and (max_map_channel[i][j]==(channel_pick+1)): 
                            
                            #if max_map_unit[i][j] not in cell_list: continue
                            #print max_map_unit[i][j], max_map_channel[i][j]
                            #if channel_pick<10: continue
                            #if channel_pick<4: continue
                            #if j != 8: continue
                            ##if depth=='cortex' and max_map_channel[i][j]<=8: continue
                            #print "Cell made thresholds: ", j

                            temp_img = max_map_list[i][j]
                            #if (max_map_nspikes[i][j]<min_spikes) and (np.max(temp_img)<0.01): continue #Skip cells with too few spikes and under 1% DF/F
                            #if np.max(temp_img) < 0.005: continue #Exclude all cells with DF/F < 1%
                            
                            #print np.max(temp_img), np.min(temp_img)
                            temp_img = ndimage.gaussian_filter(temp_img, sigma=3)
                            #temp_img = np.clip(temp_img,0,1)
                            #v_min = np.nanmin(temp_img)
                            #v_max = np.nanmax(temp_img)
                            #temp_img = (temp_img - v_min)/(v_max-v_min)
                            temp_array = np.zeros((n_pixels,n_pixels), dtype=np.float32)
                            save_indexes.append([])
                            #sum_temp0=0
                            
                            
                            #Set threhsholds for contour maps
                            ctxthresh = 0.50
                            subctxthresh = 0.50
                            delta_thresh = 0.04
                            if depth=='cortex':
                                lowthresh = np.max(temp_img)*ctxthresh
                                hightresh = np.max(temp_img)*ctxthresh+delta_thresh*np.max(temp_img)
                            else:
                                lowthresh = np.max(temp_img)*subctxthresh
                                hightresh = np.max(temp_img)*subctxthresh+delta_thresh*np.max(temp_img)

                            for k in range(len(temp_img)):
                                temp0 = np.where(np.logical_and(temp_img[k]>=lowthresh, temp_img[k]<=hightresh))[0]  #Use 1.0 for coverage maps
                                #sum_temp0+=len(temp0)
                                save_indexes[cell_counter].append(temp0)
                            #Different version of color scalling
                            colors.append((float(max_map_channel[i][j]-min(max_map_channel[i]))/(max(max_map_channel[i])-min(max_map_channel[i]))+0.2)/1.2)#16.)
                            #colors.append((channel_pick+2)/18.)
                            #colors.append(float(cell_counter+5)/10.)

                            #Overlay contour plot with jet image
                            if True:
                                plt.close()
                                temp_array = np.zeros((n_pixels,n_pixels), dtype=np.float32)
                                temp_array = np.ma.array(temp_array, mask=generic_mask_indexes, fill_value = 0, hard_mask = True)
                                temp_img = np.ma.array(temp_img, mask=generic_mask_indexes, fill_value = 0, hard_mask = True)
                                
                                #ax = plt.subplot(1,1,1)
                                #plt.imshow(temp_img, cmap=cm.jet, alpha=1)
                                #plt.title("Ch: "+str(channel_pick)+"  unit: "+str(j))
                                #plt.tick_params(
                                    #axis='both',          # changes apply to the x-axis
                                    #which='both',      # both major and minor ticks are affected
                                    #bottom='off',      # ticks along the bottom edge are off
                                    #top='off',         # ticks along the top edge are off
                                    #labelbottom='off')
                                #ax.get_xaxis().set_visible(False)
                                #ax.get_yaxis().set_visible(False)

                                channel_maps.append(temp_img)

                                #print "No. cells: ", cell_counter
                                for q in [cell_counter]:         #Compute max map from all cells in each experiment
                                    for k in range(n_pixels):   #Scan each line
                                        temp_array[k][save_indexes[q][k]]=colors[q]
                                        #temp_array[k][save_indexes[q][k]]=0.9#colors[q]
                                        #print save_indexes[k]
                                #print temp_array
                                #plt.imshow(temp_array,  cmap=cm.jet, vmin=0, vmax=1.0)
                                #plt.imshow(temp_img, cmap=cm.Greys,alpha=.6)
                                
                                #plt.savefig('/home/cat/'+str(cell_counter), format='png')#, dpi=600)
                                #plt.close()
                                
                            
                            #cell_coverage[counter-1].append(sum_temp0)
                            cell_counter+=1
                            #print cell_counter
                            
                            ###Break out of loop if only considering single cells per channel/per recording
                            #cell_break=True
                            #break
                        #if cell_break: break

                if len(channel_maps)>0: 
                    print "... plotting channel: ", channel_pick
                    plt.suptitle("Ch: "+str(channel_pick))

                    for k in range(len(channel_maps)): 
                        ax = plt.subplot(2, len(channel_maps), k+1)
                        
                        #Centred imaging
                        #v_max = np.nanmax(np.abs(channel_maps[k])); v_min = -v_max
                        #plt.imshow(channel_maps[k], vmin=v_min, vmax=v_max, cmap=cm.jet, alpha=1)

                        #v_max = np.nanmax(np.abs(channel_maps[k])); v_min = -v_max
                        plt.imshow(channel_maps[k], cmap=cm.jet, alpha=1)

                        
                        ax.get_xaxis().set_visible(False)
                        ax.get_yaxis().set_visible(False)                

                        ch_ave = np.mean(channel_maps, axis=0)
                        ax = plt.subplot(2, len(channel_maps), k+len(channel_maps)+1)


                        plt.imshow(channel_maps[k]-ch_ave, cmap=cm.jet, alpha=1)
                        ax.get_xaxis().set_visible(False)
                        ax.get_yaxis().set_visible(False)         

                    #plt.suptitle(str(v_min)+ "   " + str(v_max))
                    plt.savefig('/home/cat/'+str(channel_pick), format='png')#, dpi=600)
                    plt.show()                            
                else: 
                    print "... channel does not have map: ", channel_pick

            #for j in range(len(max_map_list[i])):       #Compute max map from all cells in each experiment
            temp_array = np.zeros((n_pixels,n_pixels), dtype=np.float32)
            temp_array = np.ma.array(temp_array, mask=generic_mask_indexes, fill_value = 0, hard_mask = True)
            #print "No. cells: ", cell_counter
            for q in range(cell_counter):         #Compute max map from all cells in each experiment
                for k in range(n_pixels):   #Scan each line
                    temp_array[k][save_indexes[q][k]]=colors[q]
                    #print save_indexes[k]
            #print temp_array
            plt.imshow(temp_array, cmap=cm.jet, vmin=0, vmax=1.0)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            plt.title(depth+" "+state+"\n #cells: "+str(cell_counter))

            #Compute average cell overlap for each state:
            if cell_counter>0:
                ave_coverage = []
                for q in range(cell_counter): #Loop over all cells - excluding MUA cell which comes last
                    cell_area = np.zeros((n_pixels,n_pixels),dtype=np.float32)
                    for i in range(n_pixels):
                        cell_area[i][save_indexes[q][i]]=1.0
                    ave_coverage.append(np.sum(cell_area)/np.sum(cortex_pixels))
                print depth, " ", state, "  coverages: ", ave_coverage
            else:
                print "Insufficient single units"

            counter+=1

            #quit()

    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    plt.suptitle(file_name, fontsize = 30)
    plt.show()


#Plot specificity contour maps
if specificity_maps_each_channel:
    ''' This routine requires loading a single SUA sort and the same file but with MUA triggers
    '''
    
    #************* LOAD MASKS *****************
    n_pixels = len(max_map_list[0][0])

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
    
    cortex_pixels = n_pixels**2-np.sum(generic_mask_indexes)

    min_ptp = 40
    if 'mua' in file_name:  #MUA should plot all mixed units regardless of amplitude
        min_ptp = 20
    min_spikes = 250
    
    ratios=[]
    
    depths = ['cortex','subcortical']
    states = ['anesthetized','awake']
    #map_list contains of list of all cells in each experiment: [recording, cell#]
    for channel_pick in np.arange(0,16,1):
        counter=1
        cell_coverage=[]
        #channel_pick = 12
        for depth in depths:
            for state in states:
                print "Processing: ", depth, state               
                ax = plt.subplot(2,2,counter)
                save_indexes=[]
                colors= []
                cell_counter=0
                cell_coverage.append([])
                for i in range(len(max_map_list)):
                    for j in range(len(max_map_list[i])):       #Compute max map from all cells in each experiment
                        if (max_map_depth[i][j]==depth) and (max_map_state[i][j]==state) and (max_map_nspikes[i][j]>min_spikes) and (max_map_ptp[i][j]>min_ptp) \
                        and (max_map_channel[i][j]==channel_pick):   
                            print "Cell made thresholds: ", j

                            #print "Processing cell: ", j
                            temp_img = max_map_list[i][j]
                            temp_img = ndimage.gaussian_filter(temp_img, sigma=3)
                            #temp_img = np.clip(temp_img,0,1)
                            v_min = np.nanmin(temp_img)
                            v_max = np.nanmax(temp_img)
                            temp_img = (temp_img - v_min)/(v_max-v_min)
                            temp_array = np.zeros((n_pixels,n_pixels),dtype=np.float32)
                            save_indexes.append([])
                            sum_temp0=0
                            sum_default=0
                            for k in range(len(temp_img)):
                                #temp0 = np.where(np.logical_and(temp_img[k]>=0.75, temp_img[k]<=0.78))[0]
                                temp0 = np.where(np.logical_and(temp_img[k]>=0.75, temp_img[k]<=1.0))[0]  #Use 1.0 for coverage maps
                                sum_temp0+=len(temp0)
                                save_indexes[cell_counter].append(temp0)
                            #colors.append(float(max_map_channel[i][j])/16.)
                            colors.append(float(cell_counter+5)/10.)

                            cell_coverage[counter-1].append(sum_temp0)
                            cell_counter+=1
                            print cell_counter
                #for j in range(len(max_map_list[i])):       #Compute max map from all cells in each experiment
                temp_array = np.zeros((n_pixels,n_pixels), dtype=np.float32)
                temp_array = np.ma.array(temp_array, mask=generic_mask_indexes, fill_value =0, hard_mask = True)
                print "No. cells: ", depth, " ", state, " = ", cell_counter
                for q in range(cell_counter):         #Compute max map from all cells in each experiment
                    for k in range(n_pixels):   #Scan each line
                        temp_array[k][save_indexes[q][k]]=colors[q]
                        #print save_indexes[k]
                #print temp_array
                plt.imshow(temp_array, cmap=cm.jet, vmin=0, vmax=1.0)
                plt.title(depth+" "+state)
                counter+=1
        
                #After finishing a state and a channel for both SUA and MUA: compute overlap between areas covered (Divide by SUA area)
                #First compute union of all cell area coverage:
                if cell_counter>1:
                    cell_area = np.zeros((n_pixels,n_pixels),dtype=np.float32)
                    for q in range(cell_counter-1): #Loop over all cells - excluding MUA cell which comes last
                        for i in range(n_pixels):
                            cell_area[i][save_indexes[q][i]]=1.0

                    print np.sum(cell_area)

                    ##If comparing SUA with MUA coverages.
                    #if False:
                        #mua_area = np.zeros((n_pixels,n_pixels),dtype=np.float32)
                        #for q in [cell_counter-1]: #Loop over all cells - excluding MUA cell which comes last
                            #for i in range(n_pixels):
                                #mua_area[i][save_indexes[q][i]]=1.0
                        #print np.sum(mua_area)
                        
                        ##print len(np.where(np.logical_and(cell_area.ravel()==mua_area.ravel()) and )[0])
                        #common_area = len(np.where(np.logical_and(cell_area==1,mua_area==1))[0])#, mua_area==1))[0]
                    
                        #ratio = common_area/np.sum(cell_area)
                        #print "overlap ratio**************************************: ", ratio
                        #ratios.append(ratio)
                else:
                    print "Insufficient single units"
                #quit()
            
            
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        plt.suptitle(file_dir, fontsize = 30)
        #plt.show()


    print ratios
    #print cell_coverage[0]
    #print "Brain covers No. pixels: ", sum_default
    #for i in range(len(cell_coverage)):
        #print np.average(cell_coverage[i])/(cortex_pixels)*100.

#Plot specificity contour maps
if specificity_maps_overlap_adjacent_channel:
    ''' This routine requires single recording loading - compares overlap of maps between cells on adjacent channels (if such pairs exist)
    '''

    #************* LOAD MASKS *****************

    plotting = True

    n_pixels = 256#len(max_map_list[0][0])

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
    
    cortex_pixels = n_pixels**2-np.sum(generic_mask_indexes)

    min_ptp = 1
    if 'mua' in file_name:  #MUA should plot all mixed units regardless of amplitude
        min_ptp = 20
    min_spikes = 1
    
    #depth = 'cortex'
    #depth = 'subcortical'
    #state = 'anesthetized'
    #state = 'awake'
    
    n_plots = 1 #Number of rows of plots
    all_overlap = []
    font_size = 20
    cell_list =[62,69]
    
    for i in range(len(max_map_list)):
        print "Processing recording: ", i
        ratios=[]
        plot_counter = 0

        overlap_areas = []
        channel_image = []
        cell_map = []
        cell_map_all = []
        ch_colors = []
        temp_imgs = []
        for channel_pick in np.arange(0,16,1):
            #print "Processing channel: ", channel_pick
            cell_map.append([])
            channel_image.append([])
            ch_colors.append([])
            temp_imgs.append([])
            cell_map_all.append([])
            
            for j in range(len(max_map_list[i])):       #Compute max map from all cells in each experiment
                cell_coverage=[]
                save_indexes=[]
                
                if (max_map_channel[i][j]==channel_pick):   
                    #print "Cell: ", j, "  ch: ", channel_pick
                    #if (j not in cell_list): continue

                    #Pick one cell from each depth - then skip to next electrode
                    #print "Cell made thresholds: ", j, " on ch: ", channel_pick

                    temp_img = max_map_list[i][j]
                    temp_img = ndimage.gaussian_filter(temp_img, sigma=3)
                    temp_imgs[channel_pick] = temp_img
                    #temp_img = np.clip(temp_img,0,1)
                    #v_min = np.nanmin(temp_img)
                    #v_max = np.nanmax(temp_img)
                    #temp_img = (temp_img - v_min)/(v_max-v_min)
                    temp_array = np.zeros((len(temp_img),len(temp_img)),dtype=np.float32)
                    save_indexes=[]
                    
                    thresh = np.max(temp_img)*0.5
                    thresh_delta = thresh*.1
                    
                    sum_temp = []
                    for k in range(len(temp_img)):
                        temp0 = np.where(np.logical_and(temp_img[k]>=thresh, temp_img[k]<=thresh+thresh_delta))[0]  #Use 1.0 for coverage maps
                        cell_map[channel_pick].append(temp0)
                        temp1 = np.where(np.logical_and(temp_img[k]>=thresh, temp_img[k]<=1.0))[0]  #Use 1.0 for coverage maps
                        sum_temp.append(temp1)

                    cell_map_all[channel_pick].append(sum_temp)
                    colors =float(channel_pick+4)/20.
                    ch_colors[channel_pick].append(colors)
               
                    #Save channel image
                    temp_array = np.zeros((n_pixels,n_pixels), dtype=np.float32)
                    temp_array = np.ma.array(temp_array, mask=generic_mask_indexes, fill_value = 0, hard_mask = True)
                    for k in range(n_pixels):   #Scan each line
                        temp_array[k][cell_map[channel_pick][k]]=colors
                    channel_image[channel_pick].append(temp_array)

        #COMPUTE OVERLAPPING AREAS for each recording
        overlap_areas=[]
        for channel_pick in np.arange(0,14,1): #loop over all channels up to end -1
            #print ""
            #print "Cells on current_channel: ", channel_pick, "  =  ", len(cell_map_all[channel_pick])
            #print "Cells on next_channel: ", channel_pick+1, "  =  ", len(cell_map_all[channel_pick+1])
            if (len(cell_map_all[channel_pick])>0) and (len(cell_map_all[channel_pick+1])>0):  #if cell on current channel and next:
                for current_cell in range(len(cell_map_all[channel_pick])):
                    for next_cell in range(len(cell_map_all[channel_pick+1])):
                        
                        common_area = []
                        #Loop over each row of the cell_map using np.intersect1d; don't divide until after all overlap area computed
                        for k in range(len(temp_img)):
                            common_area.append(np.intersect1d(cell_map_all[channel_pick][current_cell][k],cell_map_all[channel_pick+1][next_cell][k]))

                        n_pixels_ch1 = len([y for x in cell_map_all[channel_pick][current_cell] for y in x])
                        n_pixels_ch2 = len([y for x in cell_map_all[channel_pick+1][next_cell] for y in x])
                        n_pixels_common = len([y for x in common_area for y in x])

                        temp_array = np.zeros((n_pixels,n_pixels), dtype=np.float32)
                        temp_array = np.ma.array(temp_array, mask=generic_mask_indexes, fill_value = 0, hard_mask = True)
                        for k in range(n_pixels):   #Scan each line
                            temp_array[k][common_area[k]]=0.5
                        channel_overlap = temp_array

                        #Save data for arrays
                        percent_common = float(n_pixels_common)/min(n_pixels_ch1,n_pixels_ch2)*1E2 #Divide by smaller area of the two cells;
                        overlap_areas.append(percent_common)
                        all_overlap.append(percent_common)

    
        #print overlap_areas
        print "Stats for single recording:"
        print np.average(overlap_areas)
        print np.std(overlap_areas)


    print "Stats for all recordings:"
    print np.average(all_overlap)
    print np.std(all_overlap)

    print all_overlap
    quit()
    
    #After finishing a state and a channel for both SUA and MUA: compute overlap between areas covered (Divide by SUA area)
    #First compute union of all cell area coverage:
    if False:
        if cell_counter>1:
            cell_area = np.zeros((n_pixels,n_pixels),dtype=np.float32)
            for q in range(cell_counter-1): #Loop over all cells - excluding MUA cell which comes last
                for i in range(n_pixels):
                    cell_area[i][save_indexes[q][i]]=1.0

            print np.sum(cell_area)

            mua_area = np.zeros((n_pixels,n_pixels),dtype=np.float32)
            for q in [cell_counter-1]: #Loop over all cells - excluding MUA cell which comes last
                for i in range(n_pixels):
                    mua_area[i][save_indexes[q][i]]=1.0
            print np.sum(mua_area)
            
            #print len(np.where(np.logical_and(cell_area.ravel()==mua_area.ravel()) and )[0])
            common_area = len(np.where(np.logical_and(cell_area==1,mua_area==1))[0])#, mua_area==1))[0]
        
            ratio = common_area/np.sum(cell_area)
            print "overlap ratio**************************************: ", ratio
            ratios.append(ratio)
        else:
            print "Insufficient single units"
        #quit()


#Plot specificity contour maps
if specificity_maps_over_time:

    #************* LOAD MASKS *****************
    n_pixels = len(max_map_list[0][0])

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
    
    cortex_pixels = n_pixels**2-np.sum(generic_mask_indexes)

    min_ptp = 40
    if 'mua' in file_name:  #MUA should plot all mixed units regardless of amplitude
        min_ptp = 20
    min_spikes = 200
    
    depths = ['cortex','subcortical']
    states = ['anesthetized','awake']
    #map_list contains of list of all cells in each experiment: [recording, cell#]
    cell_coverage=[]
    colors= []
    color_counter=0
    cell_counter=0
    for j in range(len(max_map_list[0])):       #Loop over cells - look for identical cells
        print "Looking for cell: ", j

        for i in range(len(max_map_list)):
            sort_index = np.argsort(max_map_unit[i])
            max_map_list[i]=np.array(max_map_list[i])[sort_index]

            if (max_map_nspikes[i][j]>min_spikes) and (max_map_ptp[i][j]>min_ptp):   
                ax = plt.subplot(5,len(max_map_list),cell_counter*len(max_map_list)+i+1)
                if i==0:
                    ax.set_ylabel("Unit: "+str(j))
                if cell_counter==0:
                    plt.title(file_names[0][i]+"\n")
                print "Recording # made thresholds: ", i

                #print "Processing cell: ", j
                temp_img = max_map_list[i][j]
                temp_img = ndimage.gaussian_filter(temp_img, sigma=4)
                #temp_img = np.clip(temp_img,0,1)
                v_min = np.nanmin(temp_img)
                v_max = np.nanmax(temp_img)
                temp_img = (temp_img - v_min)/(v_max-v_min)
                temp_array = np.zeros((n_pixels,n_pixels),dtype=np.float32)
                save_indexes=[]
                for k in range(len(temp_img)):
                    temp0 = np.where(np.logical_and(temp_img[k]>=0.75, temp_img[k]<=1.))[0]  #Use 1.0 for coverage maps
                    save_indexes.append(temp0)
                colors=float(max_map_channel[i][j])/16.
                color_counter+=1

                temp_array = np.zeros((n_pixels,n_pixels), dtype=np.float32)
                temp_array = np.ma.array(temp_array, mask=generic_mask_indexes, fill_value =0, hard_mask = True)
                for k in range(n_pixels):   #Scan each line
                    temp_array[k][save_indexes[k]]=colors
                plt.imshow(temp_array, cmap=cm.jet)
                plt.title("\n"+str(max_map_nspikes[i][j]),fontsize=14)

                ax.get_xaxis().set_visible(False)
                ax.yaxis.set_major_formatter(plt.NullFormatter())
                
        cell_counter+=1
        if cell_counter == 5:
            figManager = plt.get_current_fig_manager()
            figManager.window.showMaximized()
            plt.suptitle(file_dir, fontsize = 30)
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
            #plt.tight_layout(pad=0.1, w_pad=0.9, h_pad=1.0)            
            #ax.set_xticklabels([])
            #ax.set_yticklabels([])
            ax.set_aspect('equal')
            
            plt.show()
            cell_counter = 0
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    plt.suptitle(file_dir, fontsize = 30)
    plt.show()            
    
#Plot contour maps using lfp and/or power specgrams
if specificity_maps_lfp:

    #************* LOAD MASKS *****************
    n_pixels = len(max_map_list[0][0])
    #n_pixels = 256

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

    #bands_out = ['delta', 'theta', 'alpha', 'beta', 'gamma', 'high']
    bands_out = ['all']


    lfp_static_maps = []
    channels = np.arange(0,15,1)
    #channels = [0,5,10,15]

    counter = 1
    for band in bands_out:
        temp_array = np.zeros((n_pixels,n_pixels),dtype=np.float32)
        temp_array = np.ma.array(temp_array, mask=generic_mask_indexes, fill_value =0, hard_mask = True)
        for channel in channels:
            lfp_static_map = np.load(file_dir + file_name+'/'+file_name+'_lfpmap_band_'+band+'_channel_'+str(channel).zfill(2)+'.npy')
            #lfp_static_map = np.load(file_dir + file_name+'/'+file_name+'_powermap_band_'+band+'_channel_'+str(channel).zfill(2)+'.npy')
    
            print "Processing: ", channel, " ", band             

            ax = plt.subplot(3,2,counter)
            #ax = plt.subplot(1,1,counter)
            save_indexes=[]

            #plt.imshow(lfp_static_map)
            #plt.show()
            
            colors= []
            color_counter=0
            temp_img = lfp_static_map
            temp_img = ndimage.gaussian_filter(temp_img, sigma=4)
            v_min = np.nanmin(temp_img)
            v_max = np.nanmax(temp_img)
            temp_img = (temp_img - v_min)/(v_max-v_min)
            #plt.imshow(temp_img)
            #plt.show()
            #temp_array = np.zeros((n_pixels,n_pixels),dtype=np.float32)

            for k in range(n_pixels):
                temp0 = np.where(np.logical_and(temp_img[k]>=0.75, temp_img[k]<=0.77))[0]
                save_indexes.append(temp0)
                #print temp0
            colors= float(channel)/16.
            print colors
            #plt.imshow(temp_array)
            #plt.show()


            #color_counter+=1
            #for q in range(color_counter-1):         #Compute max map from all cells in each experiment
            for k in range(n_pixels):   #Scan each line
                temp_array[k][save_indexes[k]]=colors

            plt.imshow(temp_array, cmap=cm.jet)
            plt.title(band)
        counter+=1
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    plt.suptitle(file_dir+file_name, fontsize = 30)
    plt.show()
    
if max_maps:
    ''' Average pre- and post- spike intervals to obtain static map; See also the "max" version of this function'''

    n_pixels = len(max_map_list[0][0])
    
    #Interpolate static maps to 256 x 256 pixel maps
    temp_list_max = list(max_map_list)
    temp_list_min = list(min_map_list)
    if n_pixels==128:
        for i in range(len(max_map_list)):
            for j in range(len(max_map_list[i])):
                temp_list_max[i][j]= scipy.misc.imresize(max_map_list[i][j],2., interp='bicubic', mode=None)
                temp_list_min[i][j]= scipy.misc.imresize(min_map_list[i][j],2., interp='bicubic', mode=None)
    
    max_map_list = list(temp_list_max)
    min_map_list = list(temp_list_min)
    
    n_pixels=256
    print np.array(max_map_list[0]).shape


    #Load General mask (removes background)
    generic_mask_file = []
    generic_mask_file = main_dir + 'genericmask.txt'
    if (os.path.exists(generic_mask_file)==True):
        generic_coords = np.loadtxt(generic_mask_file)


    #Load Bregma coords for computing contralateral areas
    #bregma_mask_file = file_dir + 'bregmamask.txt'
    #bregma_coords = np.loadtxt(bregma_mask_file)
    #if len(bregma_coords)==2:
        #bregma_coords = [bregma_coords]
    #bregma_coords_temp = []
    #for j in range(len(bregma_coords)):  #Remove centreline from all images
        #for k in range(n_pixels):
            #for l in range(7):
                #bregma_coords_temp.append([k,min(n_pixels-1,bregma_coords[j][1]-3+l)])

    generic_mask_indexes=np.zeros((n_pixels,n_pixels))
    for i in range(len(generic_coords)):
        generic_mask_indexes[generic_coords[i][0]][generic_coords[i][1]] = False
        
    #for i in range(len(bregma_coords_temp)):
        #generic_mask_indexes[bregma_coords_temp[i][0]][bregma_coords_temp[i][1]] = True       
        

    #Parameters set
    depths = ['cortex','subcortical']
    states = ['anesthetized','awake']
    min_ptp = 35 #minimum PTP amplitude on units 

    limit_spikes = [0, 100,500,1000,1500,2000,2500,3500,4500,7500]#,10000,15000,20000,25000]
    #limit_spikes = [3500,4500,7500]#,10000,15000,20000,25000]

    plt.figure(figsize = (8,len(limit_spikes)))
    gs1 = gridspec.GridSpec(8, len(limit_spikes))
    gs1.update(wspace=0.025, hspace=0.05) # set the spacing between axes. 

    counter=0
    plot_counter=0
    for depth in depths:
        for state in states:
            
            print "Depth: ", depth, " state: ", state
            
            #************* COMPUTE MAX MAPS *****************
            save_max_maps = []
            #Loop over list of # of spikes thresholds
            for q in limit_spikes:
                img_out = []
                
                #map_list contains of list of all cells in each experiment: [recording, cell#]
                for i in range(len(max_map_list)):
                    temp_array = []
                    #print len(max_map_list[i])
                    for j in range(len(max_map_list[i])):       #Compute max map from all cells in each experiment
                        temp_array.append(max_map_list[i][j])
                        #print max_map_nspikes[i][j],max_map_ptp[i][j],max_map_depth[i][j],max_map_state[i][j]
                    
                    temp_array = np.float32(temp_array)

                    for j in range(len(temp_array)):    #Loop over all cells in experiment and look for matching conditions
                        if (max_map_nspikes[i][j]>q) and (max_map_ptp[i][j]>=min_ptp) and (max_map_depth[i][j]==depth) and (max_map_state[i][j]==state):
                            temp=temp_array[j]
                            temp = np.ma.array(temp, mask=generic_mask_indexes, fill_value = 0.)
                            v_min2 = np.nanmin(temp)
                            v_max2 = np.nanmax(temp)
                            temp_img2 = (temp - v_min2)/(v_max2-v_min2)
                            temp_array[j]= temp_img2
                        else:
                            temp_array[j]=0.0
                    
                        #plt.imshow(temp_array[j])
                        #plt.show()
                    
                    if len(temp_array)>0:
                        temp_img2 = np.amax(temp_array,axis=0)      #matthieu method of looking for max method across all slices
                        img_out.append(temp_img2)

                #DYNAMIC NORMALIZED MAP
                img_dynamic=[]
                for i in range(len(img_out)):
                    temp_array=np.float32(img_out[i])
                    #plt.imshow(temp_array)
                    #plt.show()
                    
                    temp_array = np.ma.array(temp_array, mask=generic_mask_indexes, fill_value = 0.)
                    v_min2 = np.nanmin(temp_array)
                    v_max2 = np.nanmax(temp_array)
                    #print v_min2, v_max2
                    temp_img2 = (temp_array - v_min2)/(v_max2-v_min2)
                    img_dynamic.append(temp_img2)

                #ax=plt.subplot(8,len(limit_spikes),counter+1)
                temp_array = np.float32(img_dynamic)
                temp_array = np.amax(temp_array,axis=0)
                #plt.imshow(temp_array)
                #plt.show()

                save_max_maps.append(temp_array)
                
                ax = plt.subplot(gs1[counter])
                plt.axis('on')
                ax.set_xticklabels([])
                ax.set_yticklabels([])
                ax.set_aspect('equal')
                ax.get_xaxis().set_visible(False)
                #ax.get_yaxis().set_visible(False)
                
                if q==0: 
                    plt.ylabel(depth+'\n'+state)
                plt.imshow(temp_array)
    
                #plt.title(str(q), fontsize = 8)

                counter+=1
            
            if True:
                temp_img = save_max_maps[0]
                for i in range(1, len(save_max_maps),1):
                    temp_img = np.hstack((temp_img,save_max_maps[i]))
                
                np.save(file_dir+'max_maps_'+depth+"_"+state, temp_img)
                
                #plt.close()
                #temp_img = np.load(file_dir+'max_maps_'+depth+"_"+state+'.npy')
                #plt.imshow(temp_img)
                #plt.show()

            #************* COMPUTE MIN MAPS *****************
            save_min_maps = []
           
            #Loop over list of # of spikes thresholds
            for q in limit_spikes:
                img_out = []
                #print "Threshold: ",q
                #map_list contains of list of all cells in each experiment: [recording, cell#]
                for i in range(len(min_map_list)):          #Loop over all single cell min maps for all recordings
                    temp_array = []
                    for j in range(len(min_map_list[i])):       #Compute min map from all cells in each recordings
                        temp_array.append(min_map_list[i][j])
                    
                    temp_array = np.float32(temp_array) #Convert list of minmaps to numpy array
                    if True: #Normalize individual cell maps
                        temp_array2=[]
                        for j in range(len(temp_array)):
                            
                            if (min_map_nspikes[i][j]>=q) and (min_map_ptp[i][j]>min_ptp) and (min_map_depth[i][j]==depth) and (min_map_state[i][j]==state):
                                #print "Cell passed threshold: ", j
                                temp=np.float32(temp_array[j])
                                #temp = np.ma.array(temp, mask=generic_mask_indexes, fill_value = nan)
                                v_min2 = np.nanmin(temp)
                                v_max2 = np.nanmax(temp)
                                temp = (temp - v_min2)/(v_max2-v_min2)
                                
                                temp_array2.append(temp)
                            #else:
                            #    temp_array[j]=np.max(temp)
                    
                                #counter2+=1
                    if len(temp_array2)>0:
                        temp_array2=np.array(temp_array2)
                        temp_img2 = np.amin(temp_array2,axis=0)
                        img_out.append(temp_img2)
                    #plt.imshow(img_out[i])
                    #plt.show()
                
                #DYNAMIC NORMALIZED MAP
                img_dynamic=[]
                for i in range(len(img_out)):
                    temp_array=np.float32(img_out[i])
                    temp_array = np.ma.array(temp_array, mask=generic_mask_indexes, fill_value = nan)
                    v_min2 = np.nanmin(temp_array)
                    v_max2 = np.nanmax(temp_array)
                    temp_img2 = (temp_array - v_min2)/(v_max2-v_min2)
                    img_dynamic.append(temp_img2)

                #ax=plt.subplot(8,len(limit_spikes),counter+1)
                if len(img_dynamic)>0:
                    temp_array = np.float32(img_dynamic)
                    temp_array = np.amin(temp_array,axis=0)
                else:
                    temp_array=np.zeros((n_pixels,n_pixels),dtype=np.float32)
                
                ax = plt.subplot(gs1[counter])
                plt.axis('on')
                ax.set_xticklabels([])
                ax.set_yticklabels([])
                ax.set_aspect('equal')
                ax.get_xaxis().set_visible(False)
                #ax.get_yaxis().set_visible(False)
                #plt.title(str(q), fontsize = 8)
                plt.imshow(temp_array)

                if q==0: 
                    plt.ylabel(depth+'\n'+state)

                counter+=1

                save_min_maps.append(temp_array)

            if True:
                temp_img = save_min_maps[0]
                for i in range(1, len(save_min_maps),1):
                    temp_img = np.hstack((temp_img,save_min_maps[i]))
                
                np.save(file_dir+'min_maps_'+depth+"_"+state, temp_img)
    
    plt.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=0.01)
    plt.suptitle(file_dir, fontsize=20)
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    #plt.show()
    plt.savefig(file_dir+"total_map", format='png', dpi=600)
    plt.close()



