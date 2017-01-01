#Compute spike triggered maps from continuous ephys/gcamp recs

import sys
#sys.path.insert(0, '/home/cat/code/ephys')
from sta_map_short_utils import *

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
import matplotlib.animation as animation
import glob
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
#plt.register_cmap(cmap=blue_red1)

#**** DATA FILES
#main_dir = '/media/cat/12TB/in_vivo/tim/cat/2016_05_27_gcamp/'

main_dir = '/media/cat/12TB/in_vivo/tim/cat/2016_07_11_vsd/'

track_file = 'tr2'
img_file = '/tif_files/'+track_file+'.npy'
camera_file = '/camera_files/'+track_file+'_onoff.npy'
camera_file_pulses = '/camera_files/'+track_file+'_pulses.npy'
camera_pulse_times = '/camera_files/'+track_file+'_pulses_times.npy'

#**** PARAMETERS
overwrite = False    #overwrite previous results
n_procs = 5

window = 3          #number of secs frames pre- and post- spiking to compute STMTDs for
img_rate = 150
min_spikes = 0  #Min number of spikes in unit 
spike_mode = 'all'     #Trigger off all spkes
#spike_mode = 'lockout'     #Trigger off all spkes
#spike_mode = 'burst'   #Trigger only off 1st spike in burst, lockout following 0.5 seconds
#spike_mode = '1sec'    #Trigger only on spikes with at least 0.5sec silence on both sides

#**** LOAD SINGLE UNITS
work_dir = main_dir + 'tsf_files/'
Sort = Loadptcs(track_file, work_dir, 1, save_timestamps=False) #Auto load flag for Nick's data
Sort.filename=track_file
Sort.directory=main_dir
print "# units: ", len(Sort.units)
n_units = len(Sort.units)


#**** LOAD IMAGING AND IMAGING TIMES
if False: 
    images_raw = np.load(main_dir + img_file)
    images_raw = np.float32(images_raw.reshape(len(images_raw)/128/128, 128 ,128))
    print images_raw.shape

    #Load camera on and off times
    camera_onoff = np.load(main_dir + camera_file)
    print len(camera_onoff)

    #Find start and end of imaging in ephys data; find where camera val goes to '1'
    indexes = np.where(camera_onoff==1)[0]
    start_offset = float(indexes[0])/Sort.samplerate
    end_offset = float(indexes[-1])/Sort.samplerate

    #Load camera triggers
    camera_pulses = np.int32(np.load(main_dir + camera_file_pulses))
    print len(camera_pulses) 
    camera_pulses = camera_pulses[indexes[0]:indexes[-1]] #Look at pulses only between ON and OFF times...

    if (os.path.exists(main_dir + camera_pulse_times)==False):
        ctr = 0
        ontimes = []
        betweentimes = []
        for k in range(len(camera_pulses)-1):
            if ((camera_pulses[k] - camera_pulses[k+1])==-1): 
                ontimes.append(k)
                betweentimes.append(k-ctr)
                ctr = k
        np.save(main_dir + camera_pulse_times[:-4], ontimes)
        np.save(main_dir + camera_pulse_times[:-4]+'_betweentimes',betweentimes)

    ontimes = np.load(main_dir + camera_pulse_times)
    #betweentimes = np.load(main_dir + camera_pulse_times[:-4]+'_betweentimes.npy')
    print "Number of camera_pulses: ", len(ontimes)

    #Compute img_times based on # of images in bin
    img_times = np.linspace(0, float(end_offset-start_offset), len(images_raw))
    print "# Interpolated img_times: ", len(img_times)

    max_len = len(ontimes)
    if len(img_times)< max_len: max_len=len(img_times)
    print "len(ontimes): ", len(ontimes), "  len(img_times): ", len(img_times)

    #Compute img_rate
    img_rate = float(len(images_raw))/float(end_offset-start_offset)
    print "img_rate: ", img_rate

    for unit in range(n_units):
        print "\n\n****************"
        print "Processing Unit: ", unit
        #REMOVE TIME TO CAMERA ON FROM EPHYS TIMES; Even for triggered data, there is still ~1sec of buffered ephys data saved
        spikes = np.array(Sort.units[unit])/Sort.samplerate - start_offset 
        channel = Sort.maxchan[unit]
        Compute_spike_triggered_average(unit, channel, spikes, images_raw, window, img_times, main_dir, track_file, overwrite, img_rate, n_procs, spike_mode)

#**** LOAD GENERIC MASK
generic_mask_file = main_dir + 'genericmask.txt'
if (os.path.exists(generic_mask_file)==True):
    generic_coords = np.loadtxt(generic_mask_file)
else:
    Define_generic_mask(images_raw, main_dir)
    generic_coords = np.int32(np.loadtxt(generic_mask_file))
    
generic_mask_indexes=np.zeros((128,128))
for i in range(len(generic_coords)):
    generic_mask_indexes[generic_coords[i][0]][generic_coords[i][1]] = True


start_time = -.2  #1 sec before t=0
end_time = +.2   #3 seconds after t=0

#**** MAKE STATIC FIGS
if True:
    #select_units = [0]
    select_units = np.arange(0,n_units,1)

    plot_img = []
    for ctr, k in enumerate(select_units):
        if len(Sort.units[k])<min_spikes: continue
        print "Loading saved .npy files: ", k
        channel = Sort.maxchan[k]
        temp_name = glob.glob(main_dir + 'stm_files/img_avg_'+track_file + '_unit'+str(k).zfill(3)+"_ch"+str(channel).zfill(3)+"_all_"+str(window)+"*")[0]
        STM = np.load(temp_name)
        
        #Apply mask
        n_pixels=128
        temp_array = np.ma.array(np.zeros((len(STM),n_pixels,n_pixels),dtype=np.float32), mask=True)
        for i in range(len(STM)):
            temp_array[i] = np.ma.array(STM[i], mask=generic_mask_indexes, fill_value = 0., hard_mask = True)
        
        ax = plt.subplot(len(select_units), 1, ctr+1)
        img_out = []
        block_save = 1
        img_rate = float(len(STM))/6.
        #for i in range(0,len(STM),block_save):
        for i in range(int(img_rate*(3+start_time)),int(img_rate*(3+end_time)), block_save):
            img_out.append(np.ma.average(temp_array[i:i+block_save], axis=0))
        img_out = np.ma.hstack((img_out))
        
        v_abs = max(np.ma.max(img_out),-np.ma.min(img_out))
        plt.imshow(img_out, vmin = -v_abs, vmax=v_abs)
        #plt.imshow(img_out)

        plt.ylabel("U: " + str(k) + ", #spk: " + str(len(Sort.units[k])), rotation='horizontal',  labelpad=50, fontsize=10)
        ax.yaxis.set_ticks([])
        ax.xaxis.set_ticks([])

        if ctr==(len(select_units)-1): 
            plt.xlabel("Time from spike (sec)", fontsize=25)
            old_xlabel = np.linspace(0,img_out.shape[1], 11)
            new_xlabel = np.around(np.linspace(start_time,end_time, 11), decimals=2)
            plt.xticks(old_xlabel, new_xlabel, fontsize=18)
        
        #plt.ylabel("U: " + str(k) + ", #spk: " + str(len(Sort.units[k]))+"\nDepth: "+str(Sort.chanpos[Sort.maxchan[k]][1])+"um", fontsize=15)
        
    plt.show()
    #quit()

#**** INITIALIZE ANIMATIONS
Writer = animation.writers['ffmpeg']
writer = Writer(fps=5, metadata=dict(artist='Me'), bitrate=1800)

#select_units = [5]
select_units = np.arange(0,n_units,1)

vid_array = []
spikes_array = []
for k in select_units:
    if len(Sort.units[k])<min_spikes: continue
    print "Loading saved .npy files: ", k
    channel = Sort.maxchan[k]
    temp_name = glob.glob(main_dir + 'stm_files/img_avg_'+track_file + '_unit'+str(k).zfill(3)+"_ch"+str(channel).zfill(3)+"_"+spike_mode+"_"+str(window)+"sec*")[0]
    print temp_name
    criteria_spikes = int(temp_name[temp_name.find('spikes')-6:temp_name.find('spikes')-1])
    #if criteria_spikes < 0: continue
    print "Criteria_spikes: ", criteria_spikes
    spikes_array.append(criteria_spikes)
    STM = np.load(temp_name)

    #Apply mask
    n_pixels=128
    #temp_array = np.ma.array(np.zeros((len(STM),n_pixels,n_pixels),dtype=np.float32), mask=True)
    temp_array=[]
    for i in range(int(img_rate*(3+start_time)),int(img_rate*(3+end_time)), block_save):
    #for i in range(len(STM)):
        temp_array.append(np.ma.array(STM[i], mask=generic_mask_indexes, fill_value = 0., hard_mask = True))
    
    vid_array.append(temp_array)

#***********GENERATE ANIMATIONS
Writer = animation.writers['ffmpeg']
writer = Writer(fps=5, metadata=dict(artist='Me'), bitrate=1800)

fig = plt.figure() # make figure

print "... generating animation..." 
# make axesimage object
# the vmin and vmax here are very important to get the color map correct

#Fix size of concatenated animations
extra_row = 0
x = int(sqrt(len(vid_array)/1.4))
if x*int(x*1.4)< len(vid_array): 
    x+=1
    if (x-1)*int(x*1.4)>=len(vid_array): extra_row=-1

im=[]
for k in range(len(vid_array)):
    #ax = plt.subplot(5,8,k+1)
    ax = plt.subplot(x+extra_row,int(x*1.4),k+1)
    
    v_abs = max(np.ma.max(vid_array[k]),-np.ma.min(vid_array[k]))
    v_max = v_abs
    v_min = -v_abs
    #v_max=np.ma.max(vid_array[k])
    #v_min=np.ma.min(vid_array[k])
    #v_max=0.10
    #v_min=-0.10
    
    #print v_max, v_min
    
    ax.get_xaxis().set_visible(False)
    ax.yaxis.set_ticks([])
    ax.yaxis.labelpad = 0
    #ax.set_ylabel("",fontsize=6)
    ax.set_title(str(spikes_array[k])+", "+str(round(np.ma.min(vid_array[k])*100,2))+"%.."+str(round(np.ma.max(vid_array[k])*100,2))+"%", fontsize=5)
    
    im.append([])
    im[k] = plt.imshow(vid_array[k][0], cmap=plt.get_cmap('jet'), vmin=v_min, vmax=v_max, interpolation='none')#, vmin=0, vmax=v_max)
    #im[k] = plt.imshow(vid_array[k][0], cmap=blue_red1, vmin=v_min, vmax=v_max, interpolation='none')#, vmin=0, vmax=v_max)

#function to update figure
def updatefig(j):
    print j,
    plt.suptitle("Track: "+ track_file+ " spike mode: "+spike_mode +"\nFrame: "+str(j)+"  " +str(round((float(j)/img_rate+start_time),2))+"sec")

    # set the data in the axesimage object
    for k in range(len(vid_array)):
        im[k].set_array(vid_array[k][j])

    # return the artists set
    return im
    
# kick off the animation
ani = animation.FuncAnimation(fig, updatefig, frames=range(len(vid_array[0])), interval=100, blit=False, repeat=True)

if True:
#if save_animation:
    ani.save(main_dir+'movie_files/'+track_file+'_'+spike_mode+'_'+str(window)+'sec_units'+str(len(vid_array))+'.mp4', writer=writer)

plt.show()



quit()

#**** SAVE DATA
Save_time_course(unit, channel, spikes, Max_plot, Min_plot, Max_index, Min_index, window, len_frame, file_dir, file_name, area_names, sides, spike_mode)
