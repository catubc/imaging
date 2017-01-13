from tsf_ptcs_classes import *
from distributions import *
from utils_dongsheng import *
#from sequential_firing import *

import subprocess
from scipy import stats, signal
from scipy import interpolate
import numpy as np
import time, math
import sys
import os.path
import multiprocessing as mp
from matplotlib import animation
from matplotlib.path import Path
import glob #Wildcard searching in dir listing
import skimage
from skimage import data
from skimage.transform import rotate
            
from pylab import *
from scipy.interpolate import interp1d
import struct, array, csv
import scipy.optimize 
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline

import matplotlib.mlab as mlab
from scipy import ndimage

import PIL
from PIL import ImageFont
from PIL import Image
from PIL import ImageDraw

from libtfr import *


def Load_images(file_dir, file_name):
    
    if (os.path.exists(file_dir + file_name + '/' + file_name + '_images.npy')==False):

        print file_dir + file_name + '/' + file_name+'.tif'
        img = Image.open(file_dir + file_name + '/' + file_name+'.tif')

        counter=0
        if True:
            while True:
                try:
                    img.seek(counter)
                except EOFError:
                    break
                counter+=1
                print counter
        img.seek(0)

        #Default pic sizes
        n_pixels = img.size[0]

        #Initialize 3D image array
        n_frames = counter
        #n_frames = 69093

        images_raw = np.zeros((n_frames, 256, 256), dtype = np.float16)

        print "n_frames: ", n_frames
        for i in range(0, n_frames,1): 
            try:
                img.seek(i)
                print "Loading frame: ", i
                #images_raw [i] = np.flipud(np.fliplr(np.float16(img))) #FLIP IMAGES FOR Experiments Nov and Dec 2015
                if n_pixels==128:
                    images_raw[i] = scipy.misc.imresize(np.float16(img),2., interp='bicubic', mode=None)
                else:
                    images_raw[i] = np.float16(img) #2016-1-11 2016-1-14 experiment no flipping needed
                
                #if i>3708:# and i%10==0: 
                    #im = plt.imshow(images_raw[i])
                    #plt.title("Frame: " +str(i))
                    #figManager = plt.get_current_fig_manager()
                    #figManager.window.showMaximized()
                    #plt.show()

            except EOFError:
                break

        images_start= 0
        images_end = n_frames
        #For imaging without blank periods
        #images_start = 0
        #images_end = len(images_raw)
        
        images_raw=images_raw[images_start:images_end]

        print "Saving imaging array..."

        np.save(file_dir + file_name + '/' + file_name + '_images', images_raw)
        np.savetxt(file_dir + file_name + '/' + file_name + '_images_start_'+str(images_start)+'_end_'+
        str(images_end), [images_start, images_end])

        images_raw = np.load(file_dir + file_name + '/' + file_name + '_images.npy')
        images_raw = np.float16(images_raw)
        
        return images_raw

    else:
        images_raw = np.load(file_dir + file_name + '/' + file_name + '_images.npy')
        images_raw = np.float16(images_raw)



    if len(images_raw[0])==256:
        return images_raw

    else:
        #temp_img = np.float64(images_raw) #load stim data
        temp=[]
        for k in range(len(images_raw)):
            print k, len(images_raw)
            temp.append(scipy.misc.imresize(np.float16(images_raw[k]),2., interp='bicubic', mode=None)) #Interpolate
        
        images_raw = np.float16(temp)
        np.save(file_dir + file_name + '/' + file_name + '_images', images_raw)
        
        return np.array(temp)
                    
                    
        
def Rotate_images(images_raw, file_dir, file_name, overwrite_shifted):

    temp_name = file_dir + file_name + '/' + file_name + '_images_rotated'
    
    if overwrite_shifted or (os.path.exists(temp_name + '.npy')==False):
        #Recenter/shift images:
        print "Rotating images"
        with open(file_dir + "image_shift.txt", "r") as f: #Text file contains image shift and rotation angle info
            data = csv.reader(f)
            temp = []
            for row in data:
                temp.append(int(row[0]))
        x_shift = temp[0]
        y_shift = temp[1]
        angle = temp[2]
        
        n_pixels = len(images_raw[0])
        
        shift_img = np.float64(images_raw)

        #Rotate image
        print "Rotate angle: ", angle
        if angle != 0:

            for i in range(len(shift_img)):
                print "Rotating img: ", i
                shift_img[i] = skimage.transform.rotate(shift_img[i], angle)#, mode='constant', cval=100)
                
                #cross_line = np.zeros((n_pixels,n_pixels),dtype=np.float16)
                #cross_line[n_pixels/2,:] = 100
                #cross_line[:,n_pixels/2] = 100
                ###plt.imshow(np.hstack((shift_img[i]+cross_line,temp+cross_line)))
                #plt.imshow(shift_img[i])
                #figManager = plt.get_current_fig_manager()
                #figManager.window.showMaximized()
                #plt.show()
                ##quit()
            
            ###Plot shifted images for comparison
            #blank_line = np.zeros((n_pixels,3),dtype=np.float16)
            #cross_line = np.zeros((n_pixels,n_pixels),dtype=np.float16)
            #cross_line[n_pixels/2,:] = 100
            #cross_line[:,n_pixels/2] = 100
            #temp_plot = shift_img[5000]+cross_line
            #temp_plot2 = images_raw[5000] + cross_line
            #plt.imshow(np.hstack((temp_plot2,blank_line,temp_plot)))
            #mng = plt.get_current_fig_manager()
            #mng.resize(*mng.window.maxsize())
            #plt.title("Original data (left)    Shifted data (right)")
            #plt.show()
            #quit()

        print "...done rotating."
        np.save(temp_name,np.array(shift_img,dtype=np.float16))

    else:
        shift_img = np.load(temp_name+'.npy')

    n_pixels = len(shift_img[0])

    return shift_img
            

def Reduce_images_128pixels(images_raw):
    print images_raw.shape
    
    from skimage.measure import block_reduce
    reduced_frames = np.zeros((len(images_raw),128,128),dtype=np.float32)
    print reduced_frames.shape
    for i in range(len(images_raw)):
        print "Block reducing frame: ", i
        reduced_frames[i] = block_reduce(images_raw[i], block_size=(2,2), func=np.mean)

    return reduced_frames

def Luminance_contrast(frame_in, frame_out):
    for i in range(0,len(frame_in),1):
        for j in range(0,len(frame_in[0]),1):
            frame_out[i,j] = frame_in[i][j]/np.mean(frame_in)   #Definition of luminance-contrast from Ringach 2002
    return frame_out
    
def Baseline_regression(images_raw):        
    print "Removing baseline over all data"
    baseline = np.mean(images_temp, axis=0)
    images_temp = (images_temp - baseline)/baseline


def Load_images_start_end(file_dir, file_name, images_raw):
    
    mcd_file = file_dir + file_name+ '/' + file_name + '.mcd'
    if (os.path.exists(mcd_file)==True):
        data = MCD_read_imagingtimes_old(mcd_file)

        print "Finding beginning and end of imaging trigger on channel 17th .mcd file"
        temp_data = []
        for i in range(data['extra'].item_count):
            temp_data.append(data['extra'].get_data(i)) 
        temp_data = np.array(temp_data)[:,0]    #Select time column only from 17th channel

        start_array = []
        end_array = []
        start_array.append(temp_data[0])
        for i in range(1,len(temp_data)-1,1):
            if temp_data[i+1]-temp_data[i]>1.0:
                end_array.append(temp_data[i])
                start_array.append(temp_data[i+1])
        end_array.append(temp_data[-1])

        #Visualize different epoch starts/ends
        #print temp_data
        #for i in range(len(start_array)):
            #plt.plot([start_array[i],start_array[i]],[0,len(temp_data)])
            #plt.plot([end_array[i],end_array[i]],[0,len(temp_data)])
            #plt.axvspan(start_array[i], end_array[i], color='red', alpha=0.35)
        #plt.show()
        #quit()

    else:
        #Load epoch.txt file
        epoch_file = file_dir + file_name+ '/epochs.txt'
        data = np.loadtxt(epoch_file)
        start_array = data[:,0]
        end_array = data[:,1]
        #print start_array
        #print end_array
        

    #print start_array
    #print end_array
    #quit()

    #Load index of recording for series recordings
    rec_index_file_name = file_dir + file_name + '/rec_index.txt'
    if (os.path.exists(rec_index_file_name)==True):
        rec_index = np.loadtxt(rec_index_file_name)
        print "Recording index: ", rec_index

    #Find star/end for multi-part recordings
    img_start = start_array[int(rec_index)-1]
    img_end = end_array[int(rec_index)-1] #temp_data[len(temp_data)-1][0]
    print "img_start: ", img_start
    print "img_end: ", img_end

    #Load frame offsets from saved file (NB: RECORDING AND CORRECT IMAGING DON"T START/END @ img_times[0] and img_times[-1]
    #import shutil
    #offset_file_name = file_dir + file_name + '/' + file_name+ '_images_start'
    #temp_file_name = glob.glob(offset_file_name+'*')
    #temp_text = temp_file_name[0].replace(offset_file_name,'')
    #t0 = temp_text.index('end')
    #offset_start = int(temp_text[1:t0-1])
    #offset_end = int(temp_text[t0+4:])
    
    n_pixels = len(images_raw[0])
    n_frames = len(images_raw)

    ##Load original .tif file to count exact frames between img_start ... img_end
    #print "n_frames saved: ", n_frames
    #temp1_file = file_dir+file_name+'/'+file_name+'_n_frames_original'
    #if (os.path.exists(temp1_file)==False):
        #temp_file = file_dir+file_name+'/'+file_name+'.tif'
        #img = Image.open(temp_file)
        #n_frames_original=0
        #while True:
            #try:
                #img.seek(n_frames_original)
            #except EOFError:
                #break
            #n_frames_original+=1
        #print "saved: ", n_frames_original
        #np.savetxt(temp1_file, [n_frames_original])
    #else:
        #n_frames_original = np.loadtxt(temp1_file)
    #print "original n_frames: ", n_frames_original
    
    #Compute correct img_start and img_end
    len_frame = float(img_end - img_start) / n_frames
    img_rate = 1. / len_frame
    img_rate_file = file_dir + file_name+ '/img_rate.txt'
    np.savetxt(img_rate_file, [img_rate])

    #len_frame = float(img_end-img_start)/n_frames_original
    #img_rate = 1./len_frame
    #img_start = img_start+offset_start*len_frame
    #img_end = img_start+n_frames*len_frame
    print "shifted img_start, img_end ", img_start, img_end
    
    img_times = []
    for i in range(n_frames):
        img_times.append(img_start+i*len_frame)
    img_times = np.array(img_times)
    
    return img_start, img_end, len_frame, img_rate, n_pixels, n_frames, img_times

def Spike_averages((args)):
    global images_temp
    temp3 = args
    return images_temp[temp3]

def Spike_averages_parallel_prespike_2sec((args)):
    global images_temp
    temp4, index = args
    sum_images = images_temp[temp4[0]]
    for i in range(1, len(temp4),1):
        temp_img = images_temp[temp4[i]] - np.mean(images_temp[temp4[i][0:len(temp4[i])/3]],axis=0) #Remove average of 1st 3rd of images
        sum_images+= temp_img
    return sum_images

def Spike_averages_parallel_prespike_3sec((args)):
    global images_temp, temp_n_pixels
      
    temp3 = args

    sum_images = np.zeros((len(temp3[0]),temp_n_pixels,temp_n_pixels), dtype=np.float32)
    for i in range(0,len(temp3),1):
        baseline = np.average(images_temp[temp3[i][0:len(temp3[i])/2]], axis=0) 
        temp_img = images_temp[temp3[i]] #Remove avg of 1st half of images
        temp_frame = (temp_img - baseline)/baseline
        sum_images += temp_frame
    
    sum_images = sum_images/len(temp3)
    
    return sum_images
    

def Spike_averages_parallel_prespike_3sec_1D((args)):
    global images_temp, temp_n_pixels

    temp3 = args

    #sum_images = np.zeros((len(temp3[0]), temp_n_pixels, temp_n_pixels), dtype=np.float32)
    #vectors = np.zeros((len(temp3), 11665408), dtype=np.float32)
    #vectors = np.zeros((len(temp3), 64, 11392), dtype=np.float32) #Later recordings sampling rates have 178 frames per 6 seconds
    #vectors = np.zeros((len(temp3), 64, 11520), dtype=np.float32) #older recordings
    vectors = np.zeros((len(temp3), 64, 19072), dtype=np.float32) #50Hz recs
    
    for i in range(0, len(temp3), 1):
        baseline = np.average(images_temp[temp3[i][0:len(temp3[i])/2]], axis=0) 
        temp_img = images_temp[temp3[i]] #Remove avg of 1st half of images
        temp_frame = (temp_img - baseline)/baseline
        
        #sum_images += temp_frame
        #temp_stack = np.ma.hstack(temp_frame)
        temp_stack = np.hstack(temp_frame)
        indexes = np.isnan(temp_stack)
        temp_stack[indexes]=0

        #vectors[i] = scipy.misc.imresize(temp_stack,.25)
        vectors[i] = scipy.ndimage.interpolation.zoom(temp_stack,.25)

    
    #sum_images = sum_images/len(temp3)
    
    return vectors



def Spike_averages_parallel_globalsignalregression_1D((args)):
    global images_temp, temp_n_pixels

    temp3 = args

    vectors = []
    for i in range(len(temp3)): 
        temp_img = images_temp[temp3[i]] #select the image stack around     #temp3 already contains a stack of 178-180 frames centred on spike
        
        temp_stack = np.ma.hstack(temp_img)
        indexes_nan = np.isnan(temp_stack)
        temp_stack[indexes_nan]=0

        vectors.append(scipy.ndimage.interpolation.zoom(temp_stack,.25))
    
    return vectors




def Spike_averages_parallel((args)):
    global images_temp
    temp4, index = args
    sum_images = images_temp[temp4[0]]
    for i in range(1, len(temp4),1):
        sum_images+= images_temp[temp4[i]]
    return sum_images

def Sum_list((temp_list)):
    global temp_window, temp_img_rate, temp_n_pixels
    temp_sum = np.zeros((int(temp_window*temp_img_rate)*2, temp_n_pixels, temp_n_pixels), dtype=np.float16)
    for i in range(len(temp_list)):
        temp_sum += temp_list[i]
    return temp_sum
    
def Compute_spike_triggered_average(unit, channel, all_spikes, window, img_rate, img_times, n_pixels, images_raw, file_dir, file_name, n_procs, overwrite, stm_types, random_flag):
    '''Computes average frame values from t=-window .. +window (usually 180 to 270 frames) '''
    
    global images_temp, temp_window, temp_img_rate, temp_n_pixels
    temp_window = window
    temp_img_rate = img_rate
    temp_n_pixels = n_pixels

    print "No. of processors: ", n_procs

    #Remove spikes outside of imaging frames
    print "Total no. spikes: ", len(all_spikes)
    temp0 = np.where(np.logical_and(all_spikes>=img_times[0]+window, all_spikes<=img_times[-1]-window))[0]
    all_spikes = all_spikes[temp0]
    print "No. spikes within imaging period: ", len(all_spikes), " firing rate: ", float(len(all_spikes))/(img_times[-1]-img_times[0]), " Hz."

    #Save only spikes that are within imaging window:
    filename_spikes = file_dir+file_name+'/unit_'+str(unit).zfill(2)+ '_channel_' + str(channel).zfill(2) + '*.csv'
    file_out = glob.glob(filename_spikes)[0]
    #if os.path.exists(file_out[:-4]+"_imagingspikes.txt")==False: 
    np.savetxt(file_out[:-4]+"_imagingspikes.txt", all_spikes)       #Save only the spikes within imaging window
    #return [],[],[] 
    
    #all_spikes = True
    #plot_string = 'all'     #Trigger off all spkes
    #plot_string = '1sec'    #Trigger only on spikes with at least 0.5sec silence on both sides
    #plot_string = 'burst'   #Trigger only off 1st spike in burst, lockout following 0.5 seconds
        
    if len(all_spikes)==0:
        return images_raw, spikes, stm_type

    #stm_types = ["burst"]  # ["all", "burst", "1sec"]    #Do "all" last to ensure that the time courses are saved also
    for stm_type in stm_types: 
        spikes = all_spikes
        
        print "\n... processing stm_type: ", stm_type
        #Check to see if images already loaded and saved as .npy
        npy_file_name = file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_'+stm_type+'_'+str(window)+'sec_window_'+str(len(spikes)).zfill(5)+"_spikes"
        stack_file_name = file_dir + file_name + '/stack1D_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_'+stm_type+'_'+str(window)+'sec_window_'+str(len(spikes)).zfill(5)+"_spikes"

        #if (overwrite) or (os.path.exists(npy_file_name+'.npy')==False):
        if (overwrite) or (os.path.exists(stack_file_name+'.npy')==False):
        #if overwrite :
            images_temp = np.array(images_raw.copy(), dtype=np.float32)

            #Remove baseline - global activity regression
            #if False:
            #    print "Removing baseline over all data"
            #    baseline = np.mean(images_temp, axis=0)
            #    #baseline = np.mean(images_temp, axis=0, dtype=np.float32)
            #    images_temp = (images_temp - baseline)/baseline
            #else: print "Not removing baseline"
            
            #Use allspikes
            if stm_type =='all':
                #Use all spikes
                print "# Spikes ", len(spikes)

                temp3 = []
                for i in spikes:        #temp3 contains all the frame indexes from img_times for each spike in raster; e.g. 180 frames for each spike automatically aligned
                    temp3.append(np.where(np.logical_and(img_times>=i-window, img_times<=i+window))[0][0:2*int(window*img_rate)]) #Fixed this value or could be off +/-1 frame
            
            #Use only spikes surrounded by 500ms silence on both sides...
            elif stm_type=='1sec':
                temp5 = []
                temp5.append(spikes[0])
                counter = 0
                for i in range(1,len(spikes)-1,1):
                    if ((spikes[i]-temp5[counter])>=.5) and ((spikes[i+1]-spikes[i])>=.5):
                        temp5.append(spikes[i])
                        counter+=1
                print "# Spikes isolated in 1sec windows: ", counter
                spikes = temp5

                temp3 = []
                for i in spikes:
                    temp3.append(np.where(np.logical_and(img_times>=i-window, img_times<=i+window))[0][0:2*int(window*img_rate)])

            #Use only spikes > 500ms apart
            elif stm_type=='burst':
                temp5 = []
                temp5.append(spikes[0])
                counter = 0
                for i in range(1,len(spikes),1):
                    if ((spikes[i]-temp5[counter])>=.5):
                        temp5.append(spikes[i])
                        counter+=1
                print "# Spikes beginning of bursts: ", counter
                spikes = temp5
                
                temp3 = []
                for i in spikes:
                    temp3.append(np.where(np.logical_and(img_times>=i-window, img_times<=i+window))[0][0:2*int(window*img_rate)])
                       
            
            #IF COMPUTING RANDOM/NOISE Flag
            if random_flag:
                print "... computing all motifs..."
                temp3 = []
                clipped_img_times = img_times[180:-90]  #For 30Hz rates make sure 3sec+3sec at beginning and 3sec at end
                for i in clipped_img_times:        #temp3 contains all the frame indexes from img_times for each spike in raster; e.g. 180 frames for each spike automatically aligned
                    temp3.append(np.where(np.logical_and(img_times>=i-window, img_times<=i+window))[0][0:2*int(window*img_rate)]) #Fixed this value or could be off +/-1 frame
            
                print "...# random motifs: ", len(temp3)
                
                np.savetxt(file_out[:-4]+"_imagingspikes.txt", clipped_img_times)       #Save only the spikes within imaging window
                #quit()

            
            #***************************************************************
            #GENERATE ONLY AVERAGE MOTIFS
            if True:
                #Use allspikes
                spiking_modes = ['burst', 'last', 'tonic', 'first']
                if stm_type =='modes':      
                    mode_spikes = []
                    print file_out[:-4]+"_imagingspikes_grouped_spikes.txt"
                    with open(file_out[:-4]+"_imagingspikes_grouped_spikes.txt", 'rt') as inputfile:
                        reader = csv.reader(inputfile)
                        for row in reader:
                            mode_spikes.append(np.float32(row))
                            
                    #Load spiking modes from disk
                    for m in range(len(mode_spikes)): 
                        temp3 = []
                        for i in mode_spikes[m]:        #temp3 contains all the frame indexes from img_times for each spike in raster; e.g. 180 frames for each spike automatically aligned
                            temp3.append(np.where(np.logical_and(img_times>=i-window, img_times<=i+window))[0][0:2*int(window*img_rate)]) #Fixed this value or could be off +/-1 frame
                
                        #Initialize list to capture data and store frame sequences;  Make MP Pool for spiking mode
                        print "... computing spiking mode: ", spiking_modes[m], "  #spikes: ", len(mode_spikes[m])
                        images_triggered_temp=[]
                        temp4 = []  #Contains only first spikes from bursts in cell
                        
                        pool = mp.Pool(n_procs)
                        chunks = int(len(mode_spikes[m])/n_procs) #Break up the temp3 array into n_procs that are "chunk" long each
                        for i in range(n_procs):
                            temp4.append(temp3[i*chunks: (i+1)*chunks])

                        print "Removing average of all pre spike frames - (time: -", window, "sec .. 0sec)"
                        images_triggered_temp.extend(pool.map(Spike_averages_parallel_prespike_3sec, temp4))
                        
                        pool.close()
                        print "... done "

                        #Computing averages spikes
                        print "Summing Number of chunks: ", len(images_triggered_temp)

                        temp_images = np.zeros((int(window*img_rate)*2, n_pixels, n_pixels), dtype=np.float16)
                        for i in range(len(images_triggered_temp)):
                            temp_images += images_triggered_temp[i]
                        
                        #DIVIDE BY NUMBER OF CHUNKS; Note used to be divided by number of spikes; also residue is being thrown out...
                        images_processed = temp_images/float(len(images_triggered_temp))

                        npy_file_name = file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_'+ \
                                        stm_type+'_'+str(window)+'sec_window_'+str(len(spikes)).zfill(5)+"_spikes_"+spiking_modes[m]

                        np.save(npy_file_name, images_processed)
                                
                else: 
                    #Compute all frames based on image index;
                    print "Computing images from spike averages in parallel for window: ", window, " secs ..."
                    images_triggered_temp=[]
                    temp4 = []  #Contains only first spikes from bursts in cell

                    if len(spikes) < 30:
                        pool = mp.Pool(1)
                        temp4.append(temp3)            
                    else:
                        
                        #n_procs = 1
                        pool = mp.Pool(n_procs)
                        chunks = int(len(spikes)/n_procs) #Break up the temp3 array into n_procs that are "chunk" long each
                        for i in range(n_procs):
                            temp4.append(temp3[i*chunks: (i+1)*chunks])
                            
                        #DISCARD RESIDUE: May wish to recapture this eventually.
                        #if n_procs*chunks<len(spikes):
                        #    temp_sum = images_temp[temp3[n_procs*chunks]]
                        #    for i in range(n_procs*chunks+1,len(spikes),1):
                        #        temp_sum+=images_temp[temp3[i]]
                        #    images_triggered_temp.append(temp_sum)
                        
                    if True:
                        print "Removing average of all pre spike frames - (time: -", window, "sec .. 0sec)"
                        images_triggered_temp.extend(pool.map(Spike_averages_parallel_prespike_3sec, temp4))

                    if False:
                        indices = np.arange(len(spikes))
                        print "Removing average of 2sec - (time: -3sec .. -1sec)"
                        images_triggered_temp.extend(pool.map(Spike_averages_parallel_prespike_2sec, zip(temp4,indices)))
                    
                    pool.close()
                    print "... done "

                    #Sum over all spikes
                    print "Summing Number of chunks: ", len(images_triggered_temp)

                    temp_images = np.zeros((int(window*img_rate)*2, n_pixels, n_pixels), dtype=np.float16)
                    for i in range(len(images_triggered_temp)):
                        temp_images += images_triggered_temp[i]
                    
                    #DIVIDE BY NUMBER OF CHUNKS; Note used to be divided by number of spikes; also residue is being thrown out...
                    images_processed = temp_images/float(len(images_triggered_temp))

                    npy_file_name = file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_'+stm_type+'_'+str(window)+'sec_window_'+str(len(spikes)).zfill(5)+"_spikes"

                    np.save(npy_file_name, images_processed)


            #**********************************************************
            #SAVE ALL SINGLE SPIKE MOTIFS 
            if False:
                #SLIDING WINDOW METHOD
                if True: 
                    print "Removing average of all pre spike frames - (time: -", window, "sec .. 0sec)"
                    images_triggered_temp.extend(pool.map(Spike_averages_parallel_prespike_3sec_1D, temp4))
                    
                    print "... done... making 1D stack..."
                    stack_1D = vstack((images_triggered_temp))
                    print stack_1D.shape

                    np.save(stack_file_name, stack_1D)
                    pool.close()
        
                    return 0, 0, 0 #Dummy return variables
                    
                else:    #Global average
                    print "Removing average of all frames..."
                    baseline = np.mean(images_temp, axis=0)

                    vectors = []
                    images_temp = (images_temp - baseline)/baseline #Normalize all data to DF/F using GSR 
                    
                    images_triggered_temp.extend(pool.map(Spike_averages_parallel_globalsignalregression_1D, temp4))

                    print "... done... making 1D stack..."
                    stack_1D = vstack((images_triggered_temp))
                    print stack_1D.shape

                    np.save(stack_file_name+"_globalsignalregression", stack_1D)
                    pool.close()
        
                    return 0, 0, 0 #Dummy return variables
                  
        else: 
            print "Skipping processing of images ... loading from file"
            images_processed = np.load(npy_file_name+'.npy')

        #print "Shape mean: ", images_processed.shape
    
    return images_processed, spikes, stm_type

def Compute_stimulus_maps(main_dir,file_dir, area_names):
    '''Computes map from stimulus evoked data'''
    
    #Remove spikes outside of imaging frames
    #print "Total no. spikes: ", len(spikes)
    #temp0 = np.where(np.logical_and(spikes>=img_times[0]+window, spikes<=img_times[-1]-window))[0]
    #spikes = spikes[temp0]
    #print "No. stimuli: ", len(stimulus)

    save_map = False
    
    #************ COMPUTE ROI FROM MANUAL CROPPING *****************
    if True:
        
        global fig, cid, area_coords_left, area_coords_right
        
        stim_areas = ['hindlimb', 'forelimb', 'whisker','retrosplenial','visual'] 
        n_pixels = 256
        images_temp = Show_master_stimulus_maps(main_dir, stim_areas)

        #back_img = np.load(main_dir+'2015-7-23/2015-7-23-1/2015-7-23-1_all_maxmap_unit_00_channel_01_ptp_043_spikes_08381.npy')
        back_img = np.load(main_dir+'2015-7-23/2015-7-23-1/2015-7-23-1_all_maxmap_unit_01_channel_04_ptp_034_spikes_02566.npy')

        images_temp = images_temp + back_img*10
        
        fig, ax = plt.subplots()
        ax.imshow(images_temp)
        ax.plot([128,128],[0,256], linewidth=2, color='white')
        side = 'left'
        ax.set_title("Define Location of " + area_names[0], fontsize=30)
        #cid = fig.canvas.mpl_connect('button_press_event', define_area_manual)
        cid = fig.canvas.mpl_connect('button_press_event', define_area_circle)
        
        fig.canvas.update()
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        plt.ylim(n_pixels,0)
        plt.xlim(0,n_pixels)
        plt.show()

        #Visualize annotated areas:
        plt.imshow(images_temp)
        print area_coords_left
        temp_left = np.array(area_coords_left)
        plt.plot(temp_left[:,1],temp_left[:,0],color='white',linewidth=4)
        temp_right = np.array(area_coords_right)
        plt.plot(temp_right[:,1],temp_right[:,0],color='white',linewidth=4)
        plt.show()
        
        from matplotlib import path
        #Save left side pixels
        p = path.Path(temp_left)
        all_pts = []
        for i in range(256):
            for j in range(256):
                all_pts.append([i,j])
        pts_inside = p.contains_points(all_pts)
        
        pts_left = []
        for i in range(256):
            for j in range(256):
                if pts_inside[i*256+j]==True:
                    pts_left.append([i,j])
        
        stim_type=area_names[0]
        np.save(main_dir+stim_type+"_left", pts_left)

        #Save left side pixels
        p = path.Path(temp_right)
        all_pts = []
        for i in range(256):
            for j in range(256):
                all_pts.append([i,j])
        pts_inside = p.contains_points(all_pts)
        
        pts_right = []
        for i in range(256):
            for j in range(256):
                if pts_inside[i*256+j]==True:
                    pts_right.append([i,j])
        
        stim_type=area_names[0]
        np.save(main_dir+stim_type+"_right", pts_right)
        
        return
   
    
    
    #************ COMPUTE ROIs FROM STIMULUS ********************

    #Save annotated areas:
    file_dirs1=[]
    file_names1=[]
    if ("whisker" in area_names): 
        stim_type = 'whisker'
        file_dirs1.append('/media/cat/12TB/in_vivo/tim/dongsheng/2015-7-22/') 
        #file_names1.append(['2015-7-22-9-W1', '2015-7-22-10-w2', '2015-7-22-11-w3'])
        file_names1.append(['2015-7-22-10-w2'])

        #file_dirs1.append('/media/cat/12TB/in_vivo/tim/dongsheng/2015-7-23/')  
        #file_names1.append(['2015-7-23-5-w1','2015-7-23-6-w2'])

    if ("visual" in area_names): 
        stim_type = 'visual'
        file_dirs1.append('/media/cat/12TB/in_vivo/tim/dongsheng/2015-7-22/') 
        file_names1.append(['2015-7-22-12-v1','2015-7-22-13-V2'])
        file_dirs1.append('/media/cat/12TB/in_vivo/tim/dongsheng/2015-7-23/')  
        file_names1.append(['2015-7-23-7-v1', '2015-7-23-8-v2'])

    if ("hindlimb" in area_names): 
        stim_type = 'hindlimb'
        file_dirs1.append('/media/cat/12TB/in_vivo/tim/dongsheng/2015-7-22/') 
        file_names1.append(['2015-7-22-16-HL','2015-7-22-17-HL'])
        file_dirs1.append('/media/cat/12TB/in_vivo/tim/dongsheng/2015-7-23/')  
        file_names1.append(['2015-7-23-13-HL', '2015-7-23-14-HL'])

    if ("forelimb" in area_names): 
        stim_type = 'forelimb'
        file_dirs1.append('/media/cat/12TB/in_vivo/tim/dongsheng/2015-7-22/') 
        file_names1.append(['2015-7-22-14-FL','2015-7-22-15-FL'])
        file_dirs1.append('/media/cat/12TB/in_vivo/tim/dongsheng/2015-7-23/')  
        file_names1.append(['2015-7-23-11-FL', '2015-7-23-12-FL'])

    #if ("a" in file_names[0]) or ("A" in file_names[0]): 
    #    stim_type = 'auditory'
    #    file_dirs1.append('/media/cat/12TB/in_vivo/tim/dongsheng/2015-7-22/') 
    #    file_names1.append(['2015-7-22-16-HL', '2015-7-22-17-HL'])
    #    file_dirs1.append('/media/cat/12TB/in_vivo/tim/dongsheng/2015-7-23/')  
    #    file_names1.append(['2015-7-23-13-HL', '2015-7-23-14-HL'])


    #Load recording and master bregma_loc
    master_bregma = np.loadtxt(main_dir+'bregma.txt')
    
    #y_shift = int(-master_bregma[0]+bregma_loc[0])
    #x_shift = int(-master_bregma[1]+bregma_loc[1])
    
    n_frames=250
    
    tif_files = []
    angles = []
    x_shifts = []
    y_shifts = []
    
    for i in range(len(file_dirs1)):
        #Get rotation angle
        with open(file_dirs1[i] + "image_shift.txt", "r") as f: #Text file contains image shift and rotation angle info
            data = csv.reader(f)
            temp = []
            for row in data:
                temp.append(int(row[0]))

        #Get shift in image for experiment
        bregma_loc = np.loadtxt(file_dirs1[i]+'/bregma.txt')
        y_ = int(-master_bregma[0]+bregma_loc[0])
        x_ = int(-master_bregma[1]+bregma_loc[1])
       
        for j in range(len(file_names1[i])):
            new_files = glob.glob(file_dirs1[i] + file_names1[i][j] + '/*.tif')
            tif_files.extend(new_files)
            for k in range(len(new_files)):
                angles.append(temp[2])
                x_shifts.append(x_)
                y_shifts.append(y_)

    n_repeats = len(tif_files)
    #print n_repeats, len(angles)
    images_raw = np.zeros((n_repeats, n_frames, 256, 256), dtype = np.float32)

    all_stim_file_aligned = main_dir+stim_type+'_allstimuli_aligned_image'
    
    if (os.path.exists(all_stim_file_aligned+'.npy')==False):
        print "Number of repeats: ", n_repeats
        for i in range(n_repeats):
            print "Loading stim file: ", i
            temp_name = tif_files[i] #temp_name = file_dir + file_name + '/'+f[k]+str(i+1)+'.tif'

            img = Image.open(temp_name)

            for j in range(0, n_frames,1): 
                try:
                    img.seek(j)
                    temp_img = np.float64(img) #load stim data
                    temp_img = scipy.misc.imresize(np.float16(img),2., interp='bicubic', mode=None) #Interpolate
                    #print angles[i]
                    temp_img = skimage.transform.rotate(temp_img, angles[i]) #Rotate
                    temp_img = np.hstack((temp_img[:,x_shifts[i]:], temp_img[:,0:x_shifts[i]]))
                    images_raw[i][j] = np.vstack((temp_img[y_shifts[i]:,:], temp_img[0:y_shifts[i],:]))
                    
                except EOFError:
                    break                    

        #images_raw = np.array(images_raw)

        images_raw = np.float32(images_raw)
        data = np.average(images_raw,axis=0)
        np.save(all_stim_file_aligned, data)
    else:
        data = np.load(all_stim_file_aligned+'.npy')
        
    print data.shape

    baseline = np.average(data[0:20], axis=0)
    data = (data-baseline)/baseline

    sigma_value = 1
    #data = ndimage.gaussian_filter(data, sigma=sigma_value) 

    #Load General mask (removes background)
    generic_mask_file = []
    generic_mask_file = file_dir + 'genericmask.txt'
    
    if (os.path.exists(generic_mask_file)==True):
        generic_coords = np.loadtxt(generic_mask_file)
    generic_coords = generic_coords
     
    generic_mask_indexes=np.zeros((256,256),dtype=np.float32)
    for i in range(len(generic_coords)):
        generic_mask_indexes[generic_coords[i][0]][generic_coords[i][1]] = True
    
    temp_array = []
    for i in range(0, len(data),1):
        temp = np.ma.array(data[i], mask=generic_mask_indexes, fill_value = 0., hard_mask = True)
        temp_array.append(temp)
    data = temp_array

    v_max=0.05 #np.nanmax(data[35:])
    v_min=0 #np.nanmin(np.array(data[35:]))
    
    print v_min, v_max
    
    #Show stim vids;
    if False:
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=5, metadata=dict(artist='Me'), bitrate=1800)

        fig = plt.figure() # make figure

        # make axesimage object
        # the vmin and vmax here are very important to get the color map correct
        im = plt.imshow(data[0], cmap=plt.get_cmap('jet'), vmin=v_min, vmax=v_max, interpolation='none')#, vmin=0, vmax=v_max)
        # function to update figure
        def updatefig(j):
            # set the data in the axesimage object
            im.set_array(data[j])
            plt.title(stim_type +" (repeats: "+str(n_repeats)+")" + "\nFrame: "+ str(j))
            # return the artists set
            return im,
        # kick off the animation
        ani = animation.FuncAnimation(fig, updatefig, frames=range(len(data)), interval=125, blit=False, repeat=True)

        if False:
            ani.save(work_dir+file_name+'.mp4', writer=writer)

        plt.show()
    

    #***********************Plot contour maps
    bregma_coord = np.loadtxt(file_dir + 'bregma.txt')
    print bregma_coord

    #Find midline y_location based on bregma-lambda line; OK IF BRAIN UPRIGHT
    midline = 128 

    #temp_img = np.average(data[200:250],axis=0) #Retrosplenial from 'whisker' stimulus
    temp_img = np.average(data[32:82],axis=0)
    temp_img0 =temp_img
    sigma_offset = 3
    temp_img0 = ndimage.gaussian_filter(temp_img, sigma=sigma_offset) 

    v_min = 0
    v_max = np.nanmax(temp_img0)
    print "threshold vmax: ", v_max
    threshold = v_max*.85

    save_indexes_right = []
    save_indexes_left = []
    sum_temp0=0
    for k in range(len(temp_img)):
        temp0 = np.where(np.logical_and(temp_img0[k]>=threshold, temp_img0[k]<=v_max))[0]  #Use 1.0 for coverage maps
        save_indexes_right.append(temp0)
        temp1 = []
        for p in range(len(temp0)):
            temp_out = midline - (temp0[p]- midline)
            if temp_out <256: temp1.append(temp_out)
        save_indexes_left.append(temp1)

    temp_array = np.zeros((256,256), dtype=np.float32)
    temp_array = np.ma.array(temp_array, mask=generic_mask_indexes, fill_value =0, hard_mask = True)
    left_array = temp_array.copy()
    right_array = temp_array.copy()
    for k in range(256):   #Scan each line
        temp_array[k][save_indexes_right[k]]=0.5
        right_array[k][save_indexes_right[k]]=0.5
        
        temp_array[k][save_indexes_left[k]]=0.8
        left_array[k][save_indexes_left[k]]=0.8
        
    ax = plt.subplot(121)
    plt.imshow(temp_array, cmap=cm.jet, vmin=0, vmax=1.0)

    #Save annotated areas - need to convert to 2D vector stack first:
    vector = []
    for i in range(len(left_array)):
        for j in range(len(left_array[i])):
            if left_array[i][j]>0: vector.append([i, j])
    
    np.save(main_dir+stim_type+"_left", vector)

    vector = []
    for i in range(len(right_array)):
        for j in range(len(right_array[i])):
            if right_array[i][j]>0: vector.append([i, j])
                
    np.save(main_dir+stim_type+"_right", vector)

    #************ CONTOUR MAP ONLY
    temp_img1 = ndimage.gaussian_filter(temp_img, sigma=sigma_offset) 

    save_indexes_right = []
    save_indexes_left = []
    sum_temp0=0
    for k in range(len(temp_img)):
        temp0 = np.where(np.logical_and(temp_img1[k]>=threshold, temp_img1[k]<=threshold+0.002))[0]  #Use 1.0 for coverage maps
        save_indexes_right.append(temp0)
        temp1 = []
        for p in range(len(temp0)):
            temp_out = midline - (temp0[p]- midline)
            if temp_out <256: temp1.append(temp_out)
        save_indexes_left.append(temp1)

    temp_array = np.zeros((256,256), dtype=np.float32)
    temp_array = np.ma.array(temp_array, mask=generic_mask_indexes, fill_value =0, hard_mask = True)
    left_array = temp_array.copy()
    right_array = temp_array.copy()
    for k in range(256):   #Scan each line
        temp_array[k][save_indexes_right[k]]=0.5
        right_array[k][save_indexes_right[k]]=0.5

        temp_array[k][save_indexes_left[k]]=0.8
        left_array[k][save_indexes_left[k]]=0.8

    ax = plt.subplot(122)
    plt.imshow(temp_array, cmap=cm.jet, vmin=0, vmax=1.0)

    #Save contour map as 2d vector
    vector = []
    for i in range(len(left_array)):
        for j in range(len(left_array[i])):
            if left_array[i][j]>0: vector.append([i, j])
                
    np.save(main_dir+stim_type+"_left_contour", vector)

    vector = []
    for i in range(len(right_array)):
        for j in range(len(right_array[i])):
            if right_array[i][j]>0: vector.append([i, j])
                            
    np.save(main_dir+stim_type+"_right_contour", vector)

    plt.title(stim_type)
    plt.show()
    
    

def Show_master_stimulus_maps(main_dir, area_names):
    
    stim_types = area_names
    colors = np.linspace(0.1,1.0,8)

    #Load General mask (removes background)
    generic_mask_file = main_dir + 'genericmask.txt'
    if (os.path.exists(generic_mask_file)==True):
        generic_coords = np.loadtxt(generic_mask_file)
    generic_coords = generic_coords #Convert to 128 pixel resolution of stimulus recs
    
    generic_mask_indexes=np.zeros((256,256))
    for i in range(len(generic_coords)):
        generic_mask_indexes[generic_coords[i][0]][generic_coords[i][1]] = True
    

    #Load areas/contours
    areas_left = []
    areas_right = []
    #areas_contours_left = []
    #areas_contours_right = []
    for stim_type in stim_types:
        areas_left.append(np.load(main_dir+stim_type+"_left.npy"))
        areas_right.append(np.load(main_dir+stim_type+"_right.npy"))
        #areas_contours_left.append(np.load(main_dir+stim_type+'_left_contour.npy'))
        #areas_contours_right.append(np.load(main_dir+stim_type+'_right_contour.npy'))


    temp_array = np.zeros((256,256), dtype=np.float32)
    #temp_array_contours = np.zeros((256,256), dtype=np.float32)
    
    for k in range(len(stim_types)):
        print "stim type: ", stim_types[k]
        print len(areas_left[k])
        for p in range(len(areas_left[k])):
            temp_array[areas_left[k][p][0],areas_left[k][p][1]] = colors[k]
        
        print len(areas_right[k])
        for p in range(len(areas_right[k])):
            temp_array[areas_right[k][p][0],areas_right[k][p][1]] = colors[k]+0.1

    
    temp_array = np.ma.array(temp_array, mask=generic_mask_indexes, fill_value =0., hard_mask = True)
    
    return temp_array
    #plt.imshow(temp_array)

    #plt.plot([128,128],[0,250], color='white', linewidth=3)
    #plt.show()
    #quit()
        
        
        
def remove_bregma(event):
    global bregma_coords, images_temp, n_pix
    
    if event.inaxes is not None:
        bregma_coords.append((event.ydata, event.xdata))
        for j in range(len(bregma_coords)):
            for k in range(n_pix):
                for l in range(7):
                    images_temp[100][k][min(n_pix-1,int(bregma_coords[j][1])-3+l)]=0
        
        plt.imshow(images_temp[100])
        fig.canvas.draw()
    else:
        plt.close()
        fig.canvas.mpl_disconnect(cid)

def define_area(event):
    
    global area_coords, images_temp, n_pix, win, img_r
    
    if event.inaxes is not None:
        area_coords.append((int(event.ydata), int(event.xdata)))
        #print int(event.ydata), int(event.xdata)
        for j in range(len(area_coords)):
            for k in range(3):
                for l in range(3):
                    images_temp[100][min(n_pix-1,int(area_coords[j][0])-1+k)][min(n_pix-1,int(area_coords[j][1])-1+l)]=0

        ax.imshow(images_temp[100])
        plt.show()

    else:
        #print 'Exiting'
        plt.close()
        fig.canvas.mpl_disconnect(cid)
        
def PointsInCircum(r,n=100):
    return [(math.cos(2*pi/n*x)*r,math.sin(2*pi/n*x)*r) for x in xrange(0,n+1)]
            
def on_click(event):
    
    global coords, images_temp, ax, fig, cid
    
    n_pix = len(images_temp[0])
    
    if event.inaxes is not None:
        coords.append((event.ydata, event.xdata))
        for j in range(len(coords)):
            for k in range(3):
                for l in range(3):
                    images_temp[100][min(n_pix,int(coords[j][0])-1+k)][min(n_pix,int(coords[j][1])-1+l)]=0

        ax.imshow(images_temp[100])
        #plt.show()
        fig.canvas.draw()
                    #figManager = plt.get_current_fig_manager()
                    #figManager.window.showMaximized()
    else:
        print 'Exiting'
        plt.close()
        fig.canvas.mpl_disconnect(cid)

def bregma_point(event):
    
    global coords, images_temp, ax, fig, cid, n_pix, win
    
    if event.inaxes is not None:
        coords.append((event.ydata, event.xdata))
        for j in range(len(coords)):
            for k in range(3):
                for l in range(3):
                    images_temp[100][min(n_pix,int(coords[j][0])-1+k)][min(n_pix,int(coords[j][1])-1+l)]=0

        #ax.imshow(images_temp[100], cmap = cm.Greys_r)

        #fig.canvas.draw()
        
        plt.close()
        #fig.close()
        fig.canvas.mpl_disconnect(cid)
        return

def Define_generic_mask(images_processed, main_dir):

    global coords, images_temp, ax, fig, cid
    
    images_temp = images_processed.copy()

    fig, ax = plt.subplots()

    if (os.path.exists(main_dir + 'genericmask.txt')==False):
        coords=[]

        ax.imshow(images_processed[100])#, vmin=0.0, vmax=0.02)
        ax.set_title("Compute generic (outside the brain) mask")
        #figManager = plt.get_current_fig_manager()
        #figManager.window.showMaximized()
        cid = fig.canvas.mpl_connect('button_press_event', on_click)
        plt.show()

        #******* MASK AND DISPLAY AREAS OUTSIDE GENERAL MASK 
        #Search points outside and black them out:
        all_points = []
        for i in range(len(images_processed[0][0])):
            for j in range(len(images_processed[0][0])):
                all_points.append([i,j])

        all_points = np.array(all_points)
        vertixes = np.array(coords) 
        vertixes_path = Path(vertixes)
        
        mask = vertixes_path.contains_points(all_points)
        counter=0
        coords_save=[]
        for i in range(len(images_processed[0][0])):
            for j in range(len(images_processed[0][0])):
                if mask[counter] == False:
                    images_processed[100][i][j]=0
                    coords_save.append([i,j])
                counter+=1

        fig, ax = plt.subplots()
        ax.imshow(images_processed[100])
        plt.show()
       
        genericmask_file = main_dir + 'genericmask.txt'
        np.savetxt(genericmask_file, coords_save)

        print "Finished Making General Mask"

    #else:
    #    print "Loading saved general mask"
        
        
    if (os.path.exists(main_dir + 'bregmamask.txt')==False):
        bregma_coords = []
        print "Making Bregma mask"
        ax.imshow(images_processed[100])#, vmin=0.0, vmax=0.02)
        ax.set_title("Compute bregma mask")
        #figManager = plt.get_current_fig_manager()
        #figManager.window.showMaximized()
        cid = fig.canvas.mpl_connect('button_press_event', remove_bregma)
        plt.show()

       
        bregmamask_file = main_dir + 'bregmamask.txt'
        np.savetxt(bregmamask_file, bregma_coords)

        print "Finished Bregma Mask"

    #else:
    #    print "Loading saved bregma mask"
        
    return generic_coords
    
def mouse_coords(event):
    
    global xycoords, fig, cid
    
    if event.inaxes is not None:
        print "Mouse coords: ", event.ydata, event.xdata
        xycoords = [event.ydata, event.xdata]
        plt.close()
        fig.canvas.mpl_disconnect(cid)
        
        #return coords
    #else:
    #    print 'Exiting'

def define_area_circle(event):
    
    global area_coords, fig, cid, circle_size
    
    if event.inaxes is not None:
        #Define N points on a circle centred at mouse click; shift circle to location
        points = np.vstack(PointsInCircum(circle_size,n=20))
        points = points + [int(event.ydata), int(event.xdata)]

        area_coords = []
        for i in range(len(points)):
            area_coords.append((points[i][0], points[i][1]))
            
        #Plot recent area
        temp_coords = np.array(area_coords)
        plt.plot(temp_coords[:,1],temp_coords[:,0],color='white',linewidth=4)
        
        fig.canvas.draw()
        
    else:
        plt.close()
        fig.canvas.mpl_disconnect(cid)

def define_area_rectangle(event):

    global area_coords, fig, cid, circle_size
    
    if event.inaxes is not None:
        #Define N points on a circle centred at mouse click; shift circle to location
        points = np.vstack(([0,0],[60,0], [60,12],[0,12], [0,0]))
        points = points + [int(event.ydata), int(event.xdata)]

        area_coords = []
        for i in range(len(points)):
            area_coords.append((points[i][0], points[i][1]))
            
        #Plot recent area
        temp_coords = np.array(area_coords)
        plt.plot(temp_coords[:,1],temp_coords[:,0],color='white',linewidth=4)
        
        fig.canvas.draw()
        
    else:
        plt.close()
        fig.canvas.mpl_disconnect(cid)


def define_area_circle_old(event):
    
    global area_coords, fig, cid, area_coords_left, area_coords_right
    
    if event.inaxes is not None:
        #Define N points on a circle centred at mouse click; shift circle to location
        points = np.vstack(PointsInCircum(3*(256/100),n=20))
        points = points + [int(event.ydata), int(event.xdata)]

        area_coords_left = []
        area_coords_right = []
        for i in range(len(points)):
            area_coords_left.append((points[i][0], points[i][1]))
            area_coords_right.append((points[i][0], 256-points[i][1]))
            
        #Plot recent area
        temp_coords = np.array(area_coords_left)
        plt.plot(temp_coords[:,1],temp_coords[:,0],color='white',linewidth=4)
        
        temp_coords = np.array(area_coords_right)
        plt.plot(temp_coords[:,1],temp_coords[:,0],color='white',linewidth=4)

        fig.canvas.draw()
        
    else:
        plt.close()
        fig.canvas.mpl_disconnect(cid)

def define_area_manual(event):
    
    global area_coords, images_temp, n_pix, win, img_r,fig
    
    if event.inaxes is not None:
        area_coords.append((int(event.ydata), int(event.xdata)))
        #print int(event.ydata), int(event.xdata)
        for j in range(len(area_coords)):
            for k in range(3):
                for l in range(3):
                    images_temp[int(area_coords[j][0])-1+k][int(area_coords[j][1])-1+l]=0
        
        ax.imshow(images_temp)
        fig.canvas.draw()

    else:
        plt.close()
        fig.canvas.mpl_disconnect(cid)
        
def Load_max_map(file_dir, area_name):
    
    global xycoords, fig, cid, n_pix
    n_pix = 256
    xycoords = []

    #Load normalized maps and stack them for display in proper order
    depths = ['cortex','subcortical']
    states = ['anesthetized', 'awake']
    maptypes = ['max', 'min']
    counter = 0
    for depth in depths:
        for state in states:
            for maptype in maptypes:
                temp = np.load(file_dir+maptype+'_maps_'+depth+'_'+state+'.npy')
                if counter>0:
                    img = np.vstack((img,temp))
                else:
                    img=temp
                counter+=1
    
    plt.close()
    #Display all maps and use mouseclick to select particular map
    fig, ax = plt.subplots()
    ax.imshow(img)
    ax.set_title(file_dir + "\nSelect map to search for: " + area_name, fontsize = 30)
    cid = fig.canvas.mpl_connect('button_press_event', mouse_coords)
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    plt.show()

    #print "Click coords: ", xycoords
    #print "Size of image: ", 8, len(img[0])/n_pix
    #Select correct map from mouse click
    height = 8
    width = len(img[0])/n_pix
    for h in range(height):
        if (xycoords[0]>(h*n_pix)) and (xycoords[0]<((h+1)*n_pix)):
            y_panel = h
            break
    for w in range(width):
        if (xycoords[1]>(w*n_pix)) and (xycoords[1]<((w+1)*n_pix)):
            x_panel = w
            break

    counter = 0
    exit2=False
    for depth in depths:
        for state in states:
            for maptype in maptypes:
                if (counter == y_panel): 
                    exit2 = True
                    break
                counter+=1
            if exit2: break
        if exit2: break

    #Load correct row map and correct column
    images_temp = np.load(file_dir+maptype+'_maps_' + depth+'_'+state+'.npy')  
    max_maps = []
    for i in range(len(images_temp[0])/len(images_temp)):
        max_maps.append(images_temp[:,i*n_pix:(i+1)*n_pix])
    images_temp = max_maps[x_panel]
    
    return images_temp #Return only selected map for annotation

def Define_cortical_areas(file_dir, file_name, area_names, sides):
    print "Defining cortical areas"

    #Define maps from max-min super amax maps
    if True:
        global area_coords, ax, fig, cid, circle_size #Not sure need all these vars
        n_pixels = 256
      
        depth = 'cortex'#, 'subcortical']
        
        circle_sizes = [10,10,15,15,15,15,10,8]
        
        #for depth in depths:
        counter=0
        for area in area_names:
            circle_size = circle_sizes[counter]
            
            for side in sides:
                save_file = file_dir + depth+"_"+area+"_"+ side
                if (os.path.exists(save_file+'.npy')==False):
                    
                    #print "Which map to load for marking hindlimb?"
                    
                    images_temp = Load_max_map(file_dir, depth+ " " + area+" " + side)

                    area_coords = []
                    fig, ax = plt.subplots()
                    ax.imshow(images_temp)
                    ax.set_title(file_dir+"\nDefine Location of "+depth+" " + area+" "+side+' in ', fontsize=30)
                    #cid = fig.canvas.mpl_connect('button_press_event', define_area_manual)
                    
                    if counter==3:
                        cid = fig.canvas.mpl_connect('button_press_event', define_area_rectangle)
                    else:
                        cid = fig.canvas.mpl_connect('button_press_event', define_area_circle)
                    
                    #fig.canvas.update()
                    figManager = plt.get_current_fig_manager()
                    figManager.window.showMaximized()
                    plt.ylim(n_pixels,0)
                    plt.xlim(0,n_pixels)
                    plt.show()

                    #Convert coords into x and y data for plotting; select each column, then tack on first coordinate to complete circle
                    area_coords.append(area_coords[0])
                    area_coords = np.array(area_coords)
                    
                    from matplotlib import path
                    #Compute pixel locations inside cropped area
                    p = path.Path(area_coords)
                    all_pts = []
                    for i in range(256):
                        for j in range(256):
                            all_pts.append([i,j])
                    pts_inside = p.contains_points(all_pts)
                    
                    #Generate mask for saving 
                    mask_save = np.zeros((256,256),dtype=int8)+1
                    for i in range(256):
                        for j in range(256):
                            if pts_inside[i*256+j]==True:
                                mask_save[i,j]=False
                    
                    #Save mask
                    np.save(save_file, mask_save)
                    np.save(file_dir + "subcortical_"+area+"_"+ side, mask_save)
                    #Save contour
                    np.save(save_file+'_contour', area_coords)
                    np.save(file_dir + "subcortical_"+area+"_"+ side+'_contour', mask_save, area_coords)
                    
                    
            counter+=1

    else:
        print "... not defining areas.. loading stim based areas"
    
    
       
    
def Remove_global_baseline(images_raw):
    '''Global activity regression - removes baseline average over entire recording'''

    images_temp = np.array(images_raw.copy(), dtype=np.float32)
    print "Removing baseline over all data"
    baseline = np.mean(images_temp, axis=0)
    images_out = (images_temp - baseline)/baseline

    return images_out
        
def Compute_luminance_contrast(images_baseline_removed):

    print "Computing luminance contrast"

    n_frames = len(images_baseline_removed)
    width = len(images_baseline_removed[0])
    height = len(images_baseline_removed[0][0])
    
    images_contrast = np.zeros((n_frames,width,height),dtype=np.float32)
    for f in range(n_frames):
        print "Frame: ", f
        for i in range(0,width,1):
            for j in range(0,height,1):
                images_contrast[f][i,j] = images_baseline_removed[f][i,j]/np.mean(images_baseline_removed[f])   #Definition of luminance-contrast from Ringach 2002
        #ax=plt.subplot(2,1,1)
        #plt.imshow(images_baseline_removed[f])
        #ax=plt.subplot(2,1,2)
        #plt.imshow(images_contrast[f])
        #plt.show()
        
    return images_contrast

def Parallel_recleastsq((args)):
    ''' Computes '''
    global lum_contrast_movie, movie_frames, frame_indexes
    
    tau = args
    print tau

    n_frames = len(lum_contrast_movie)
    delta = 1E6
    P = delta * (np.zeros(len(movie_frames[0]), dtype=np.float32)+1.0)
    w = np.zeros(len(movie_frames[0]), dtype=np.float32)
    N = len(movie_frames[0])*len(movie_frames[0])
    mu = 0 
    beta = 0.99
    
    frame_spikes = []       #Keep track of the number of spikes fired in each frame (+ tau)
    frame_out = np.zeros((len(movie_frames[0]), len(movie_frames[0][0])),dtype=np.float32) #Initialize once and zero every frame
    for i in range(0,n_frames,1): #Loop over movie frames
        #print "Procesing frame : ", i
        I_ijn = lum_contrast_movie[i]                                                       #Line 5 
        #I_ijn = movie.frames[i]
        
        u_n = I_ijn                                                                         #Line 6: no need to subsample unless using entire frame to start
        #compute number of spikes that fall in the frame
        r_n = len(np.where(np.logical_and(frame_indexes>=i+tau, frame_indexes <= i+tau+1))[0])      #Line 7
        #print "No. spikes: ", r_n
        frame_spikes.append(r_n)
        mu = beta*mu + (1-beta)*r_n                                                         #Line 8
        r_n = r_n - mu                                                                      #Line 9
        I_n = u_n * P                                                                       #Line 10
        kappa_n = 1 + I_n*u_n                                                               #Line 11
        k_n = (P*u_n)/kappa_n                                                               #Line 12
        
        #alpha_n = r_n -w * u_n                                                              #Line 13
        #w = w+ k_n * alpha_n                                                                #Line 14
        
        T = k_n * I_n                                                                       #Line 14
        P = P - T
    
    #Decorrelate the stimulus
    decorrelated_frames = []
    for i in range(0,n_frames,1): #Loop over movie frames
        temp = P * movie_frames[i]
        decorrelated_frames.append(temp)
    
    strf = np.zeros((len(movie_frames[0]), len(movie_frames[0][0])),dtype=np.float32)
    for i in range(0,n_frames,1):
        strf+= decorrelated_frames[i] * frame_spikes[i]
    
    return strf

def Parallel_find_previous_frame_index((targets)):
    global din_times
    indexes = []
    for i in range(len(targets)):
        index = np.argmin(np.abs(din_times - targets[i]))
        
        if (targets[i] < din_times[index]):
            indexes.append(index-1)
        else:    
            indexes.append(index)
    return indexes

def Compute_decorrelated_images(images_areas_nobaseline,file_dir, file_name, unit,channel,plot_string,window,spikes_in_window, n_procs,img_times, generic_mask_indexes):
    '''Recursive least squares on imaging data'''
    
    #global lum_contrast_movie, movie_frames
    
    global lum_contrast_movie, movie_frames, frame_indexes, din_times
       
    lum_contrast_movie = images_areas_nobaseline
    movie_frames = images_areas_nobaseline
   
    tau = np.arange(-90,90,1)  #No. of frames to look at before and after spike

    din_times = img_times

    spike_times = spikes_in_window    
    
    reclsq_file = file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_'+plot_string+'_'+str(window)+'sec_window_decorrelated'
    
    if (os.path.exists(reclsq_file+".npy")==False):
        print "Computing reclsq"
        min_spikes = 30
        frame_indexes=[]
        if len(spike_times)>min_spikes:
            ##parallel computation frame_indexes
            chunks = int(len(spike_times)/n_procs) #Break up the array into n_procs that are "chunk" long each
            temp4=[]
            for i in range(n_procs):
                temp4.append(spike_times[i*chunks:(i+1)*chunks])
            
            print "Finding previous frame indexes"
            pool = mp.Pool(n_procs)
            frame_indexes.extend(pool.map(Parallel_find_previous_frame_index, temp4))
            pool.close()
            
            frame_indexes = np.sort(np.array(frame_indexes).ravel(),axis=0)
        
        #Need to provide the 3 global variables above + tau frame index as argument
        print "Decorrelating images"
        images_decorrelated = []
        pool = mp.Pool(n_procs)
        images_decorrelated.extend(pool.map(Parallel_recleastsq,tau))
        pool.close()
        
        np.save(reclsq_file,images_decorrelated)
    else:
        print "Skipping computation of reclsq"
        images_decorrelated = np.load(reclsq_file+".npy")
    
    return images_decorrelated, tau


def remove_artifact(event):
    
    global coords, images_temp, n_pix, win, img_r
    
    if event.inaxes is not None:
        coords.append((event.ydata, event.xdata))
        for j in range(len(coords)):
            for k in range(7):
                for l in range(7):
                    xx = min(n_pix-1,int(coords[j][0])-1+k)
                    yy = min(n_pix-1,int(coords[j][1])-1+l)
                    images_temp[100][xx][yy]=0
                    if (np.array(coords) == [xx,yy]).all(-1).any(): #Check to see if pair already in list
                        pass
                    else:
                        coords.append([xx,yy])

        ax.imshow(images_temp[100])
        ax.set_title("Remove Artifacts")
        plt.show()

    else:
        print 'Exiting'
        plt.close()
        fig.canvas.mpl_disconnect(cid)

def load_tif(work_dir, file_name):
    
    if (os.path.exists(work_dir + file_name +'.npy')==False):
        print "Opening: ", work_dir + file_name+'.tif'
        img = Image.open(work_dir + file_name+'.tif')

        counter=0
        if True:
            while True:
                try:
                    img.seek(counter)
                except EOFError:
                    break
                counter+=1

        n_pixels = 128
        images_raw = np.zeros((counter, n_pixels, n_pixels), dtype = np.float32)

        for i in range(0, counter,1): 
            try:
                img.seek(i)
                images_raw [i] = img 
            except EOFError:
                break

        print "Saving imaging array..."

        #n_frames = min(counter,n_frames) #make sure you don't go past end of data; NB: technically it should be n_frames-offset_frames
        #images_raw = images_raw[offset_frames:offset_frames+n_frames]

        np.save(work_dir + file_name, images_raw)
        
    else:
        print "Loading .npy file from disk"
        images_raw = np.load(work_dir+file_name+'.npy')

    return images_raw

def remove_bregma2(event):
    global bregma_coords, images_temp, n_pixels
    
    if event.inaxes is not None:
        bregma_coords.append((event.ydata, event.xdata))
        for j in range(len(bregma_coords)):
            for k in range(n_pix):
                for l in range(7):
                    images_temp[200][k][min(n_pix-1,int(bregma_coords[j][1])-3+l)]=0
        
        plt.imshow(images_temp[200])
        fig.canvas.draw()
    else:
        plt.close()
        fig.canvas.mpl_disconnect(cid)


def remove_artifact2(event):
    global coords, images_temp, n_pix
    
    if event.inaxes is not None:
        coords.append((event.ydata, event.xdata))
        for j in range(len(coords)):
            for k in range(7):
                for l in range(7):
                    xx = min(n_pix-1,int(coords[j][0])-1+k)
                    yy = min(n_pix-1,int(coords[j][1])-1+l)
                    images_temp[200][xx][yy]=0
                    #if (np.array(coords) == [xx,yy]).all(-1).any(): #Check to see if pair already in list
                    #    pass
                    #else:
                    #    coords.append([xx,yy])

        plt.imshow(images_temp[200])
        fig.canvas.draw()

    else:
        print 'Exiting'
        plt.close()
        fig.canvas.mpl_disconnect(cid)


#def mask_data(data, work_dir, file_name):
    #print "Masking data..."
    
    #global coords, ax, n_pix, images_temp
    
    #images_temp = data

    #img_rate = 150      #Frame rate of imaging
    #n_pixels = 128      #number of pixels
    #n_pix = n_pixels

    ##DEFINE GENERAL MASK******************
    #generic_mask_file = []
    #generic_mask_file = work_dir + 'genericmask.txt'
    #if (os.path.exists(generic_mask_file)==False):
        #generic_coords = Define_generic_mask(data, work_dir, file_name)
    #else:
        #generic_coords = np.loadtxt(generic_mask_file)
    
    ##DEFINE BREGMA************************
    #if (os.path.exists(work_dir + 'bregmamask.txt')==False):
        #print "Generating bregma mask..."
        #plt.close()
        #fig, ax = plt.subplots()

        #bregma_coords = []
        #print "Making Bregma mask"
        #ax.imshow(images_raw[200])
        #ax.set_title("Compute bregma mask")
        #cid = fig.canvas.mpl_connect('button_press_event', remove_bregma2)
        #plt.show()
       
        #bregmamask_file = work_dir + 'bregmamask.txt'
        #np.savetxt(bregmamask_file, bregma_coords)
        
        #bregma_coords_temp = []
        #for j in range(len(bregma_coords)):  #Remove centreline from all images
            #for k in range(n_pixels):
                #for l in range(7):
                    #bregma_coords_temp.append([k,min(n_pixels-1,bregma_coords[j][1]-3+l)])
        #bregma_coords = bregma_coords_temp
                    
        #print "... done bregma mask..."

    #else:
        #print "Loading saved bregma mask"
        ##Load Bregma coords for computing contralateral areas
        #bregma_mask_file = work_dir + 'bregmamask.txt'
        #bregma_coords = np.loadtxt(bregma_mask_file)
        #if len(bregma_coords)==2:
            #bregma_coords = [bregma_coords]
            
        #bregma_coords_temp = []
        #for j in range(len(bregma_coords)):  #Remove centreline from all images
            #for k in range(n_pixels):
                #for l in range(7):
                    #bregma_coords_temp.append([k,min(n_pixels-1,bregma_coords[j][1]-3+l)])
        #bregma_coords = bregma_coords_temp

    #generic_mask_indexes=np.zeros((n_pixels,n_pixels))
    #for i in range(len(generic_coords)):
        #generic_mask_indexes[generic_coords[i][0]][generic_coords[i][1]] = True
    #for i in range(len(bregma_coords)):
        #generic_mask_indexes[bregma_coords[i][0]][bregma_coords[i][1]] = True       

    ##MASK THE DATA; NB: PROBABLY FASTER METHOD **************************
    #temp_array = []
    #for i in range(0, len(data),1):
        #temp = np.ma.array(data[i], mask=generic_mask_indexes, fill_value = 0., hard_mask = True)
        #temp_array.append(temp)

    #return temp_array

def load_pullfile(pull_file, images_raw, trigger_value):

    import re
    text_file = open(pull_file+'.txt', "r")

    lines = text_file.read().splitlines()
    file_text = []
    for line in lines:
        file_text.append(re.split(r'\t+',line))

    file_text = np.array(file_text)

    events = file_text[1:,1]
    pull_times = [float(x) for x in file_text[1:,3]]

    print "Computing baseline..."
    print "N frames: ", len(images_raw)
    img_rate = len(images_raw)/pull_times[-1]
    print "img rate: ", img_rate
    window = int(2 * img_rate) #Number of seconds pre and post level pull

    frame_triggers = []
    for i in range(len(events)):
        if events[i]==trigger_value:
            frame_triggers.append(int(img_rate*pull_times[i]))
    
    print "Number of pulls: ", len(frame_triggers)
    
    return pull_times, frame_triggers, window, img_rate

def filter_data():
    lowcut = 1
    highcut=70
    img_rate = 150        #Frame rate of imaging

    img_rate = len(images_raw)/pull_times[-1]
    window = int(2 * img_rate) #Number of seconds pre and post level pull
    
    data_array = np.zeros((window*2,128,128), dtype=np.float32)
    print "Computing pull triggered vids..."
    for trigger in frame_triggers:
        print "pulltime frame#: ", trigger
        data_chunk = images_raw[trigger-window:trigger+window]
        
        data_chunk = np.array(butter_bandpass_filter(data_chunk, lowcut, highcut, img_rate, order = 2))
        
        data_array+=data_chunk
    
    data = data_array/len(frame_triggers)
    


def remove_3sec_baseline(images_raw, frame_triggers, pull_times, work_dir, file_name):
    
    temp_fname = work_dir+file_name+"_3sec"
    if (os.path.exists(temp_fname+'.npy')==False):

        img_rate = len(images_raw)/pull_times[-1]
        window = int(2 * img_rate) #Number of seconds pre and post level pull
        
        data_array = np.zeros((window*2,128,128), dtype=np.float32)
        
        print "Computing pull triggered vids..."
        for trigger in frame_triggers:
            print "pulltime frame#: ", trigger
            if trigger <window: continue
            data_chunk = images_raw[trigger-window:trigger+window]
            baseline = np.average(data_chunk[0:window], axis=0)
            data_chunk = (data_chunk-baseline)/baseline
            
            #sigma_value = 1
            #data_array += ndimage.gaussian_filter(data_chunk, sigma=sigma_value) 
        
            data_array += data_chunk
            
        data = data_array/len(frame_triggers)
        
        #plt.imshow(data[0])
        #plt.show()
        
        np.save(temp_fname, data)
    
    else:
        data = np.load(temp_fname+'.npy')
        
    return data

def Define_bregma_lambda(images_rotated, main_dir, file_dir, file_name):
    
    global coords, images_temp, ax, fig, cid, n_pix, bregma_coords
    
    images_temp = images_rotated.copy()
    
    plt.close()
    n_pixels = 256
    n_pix = n_pixels
    
    if (os.path.exists(file_dir + 'bregma.txt')==False):

        coords = []
        fig, ax = plt.subplots()

        ax.imshow(images_temp[100], cmap = cm.Greys_r)
        ax.set_title("Define Bregma")
        plt.ylim(255,0)
        plt.xlim(0,255)
        for k in range(n_pixels/10):
            plt.plot([k*10,k*10],[1,n_pixels-2], color='red')
            plt.plot([1,n_pixels-2],[k*10,k*10], color='red')
        cid = fig.canvas.mpl_connect('button_press_event', bregma_point)
        plt.show()

        #Search points and black them out:
        images_temp[100][coords[0][0]][coords[0][1]]=1.0

        fig, ax = plt.subplots()
        ax.imshow(images_temp[100], cmap = cm.Greys_r)
        plt.show()
        
        bregma_file = file_dir + 'bregma.txt'
        np.savetxt(bregma_file, coords)
        #if '7-22' in file_dir:  np.savetxt(main_dir, coords)

        
    #else:
    #    print "Loading saved bregma"
    
    #if False:
        #if (os.path.exists(file_dir + 'lambda.txt')==False):

            #coords = []
            #fig, ax = plt.subplots()

            #ax.imshow(images_temp[100], cmap = cm.Greys_r)
            #ax.set_title("Define Lambda")
            #for k in range(n_pixels/10):
                #plt.plot([k*10,k*10],[0,n_pixels-1], color='red')
                #plt.plot([0,n_pixels-1],[k*10,k*10], color='red')
            #cid = fig.canvas.mpl_connect('button_press_event', bregma_point)
            #plt.show()

            ##Search points and black them out:
            #images_temp[100][coords[0][0]][coords[0][1]]=1.0

            #fig, ax = plt.subplots()
            #ax.imshow(images_temp[100], cmap = cm.Greys_r)
            #plt.show()
            
            #lambda_file = file_dir + 'lambda.txt'
            #np.savetxt(lambda_file, coords)
            
        #else:
            #print "Loading saved lambda"
            

    file_name = file_dir + file_name+'/'+file_name+'_images_aligned.npy'
    if (os.path.exists(file_name)==False):
        bregma_loc = np.loadtxt(file_dir+'/bregma.txt')
        #lambda_loc = np.loadtxt(file_dir+'/lambda.txt')
        print "exp bregma: ", bregma_loc 
        
        #Realign images:
        master_bregma = np.loadtxt(main_dir+'bregma.txt')
        print "master bregma: ", master_bregma
        
        y_shift = int(-master_bregma[0]+bregma_loc[0])
        x_shift = int(-master_bregma[1]+bregma_loc[1])
        
        #Shift image - x and y directions
        for k in range(len(images_temp)):
            temp_array = np.hstack((images_temp[k][:,x_shift:n_pixels], images_temp[k][:,0:x_shift]))
            images_temp[k] = np.vstack((temp_array[y_shift:n_pixels,:], temp_array[0:y_shift,:]))
        
        images_temp = np.float16(images_temp)
        np.save(file_name, images_temp)
        
        #plt.imshow(images_temp[100])
        #plt.title("Realigned images")
        #plt.show()
        
    else:
        images_temp = np.load(file_name)
    
        #plt.imshow(images_temp[100])
        #plt.title("Realigned images")
        #plt.show()
    
    return images_temp
   
    
    
    
def Define_artifact(images_processed, file_dir, file_name, window, n_pixels,img_rate):
    ''' Tool to manually ablate small regions of imaging area that may be 
        artifacts '''
    
    global coords, images_temp, ax, fig, cid, n_pix, win
    n_pix = n_pixels
    win = window

    print "Manual Mask Mode"
    images_temp = np.array(images_processed).copy()
    
    #Load Generic Mask
    if (os.path.exists(file_dir + 'genericmask.txt')==True):
        print "Loading existing generic mask"
        generic_mask_file = file_dir +'genericmask.txt'
        coords_generic = np.loadtxt(generic_mask_file)
        #Ablate generic map
        for i in range(len(images_temp)):
            for j in range(len(coords_generic)):
                images_temp[i][min(n_pixels-1,int(coords_generic[j][0]))][min(n_pixels-1,int(coords_generic[j][1]))]=0
                
    bregma_mask_file = file_dir + 'bregmamask.txt'
    if (os.path.exists(bregma_mask_file)==True):
        bregma_coords = np.loadtxt(bregma_mask_file)
        print "Loading bregma mask"
        #Remove centreline artifacts
        for i in range(len(images_temp)):
            for j in range(len(bregma_coords)):
                for k in range(n_pix):
                    for l in range(7):
                        images_temp[i][k][min(n_pix-1,int(bregma_coords[1])-3+l)]=0  #DON"T HARDWIRE MID SLICE

    #Load existing artiact  mask file
    coords=[]
    specific_mask_file = file_dir +'artifactmask.txt'
    if (os.path.exists(specific_mask_file)==True):
        temp_data= np.loadtxt(specific_mask_file)
        for i in range(len(temp_data)):
            coords.append(temp_data[i])
        #update_length=len(coords)

        #Ablate specific map
        for i in range(len(images_temp)):
            for j in range(len(coords)):
                for k in range(7):
                    for l in range(7):
                        images_temp[i][min(n_pixels-1,int(coords[j][0])-3+k)][min(n_pixels-1,int(coords[j][1])-3+l)]=0
    #else:
        #update_length=0
        
        fig, ax = plt.subplots()
        ax.imshow(images_temp[100])
        ax.set_title("Remove Artifacts")
        cid = fig.canvas.mpl_connect('button_press_event', remove_artifact)
        plt.show()
    
    #Save total map containing specific coords
    if len(coords)>0: np.savetxt(specific_mask_file, np.array(coords))

def Ablate_outside_area(n_pixels, contour_coords):
    ''' Function takes points from contour and returns list of coordinates
        lying outside of the contour - to be used to ablate image    '''
    
    #Search points outside and black them out:
    all_points = []
    for i in range(n_pixels):
        for j in range(n_pixels):
            all_points.append([i,j])

    all_points = np.array(all_points)
    vertixes = np.array(contour_coords) 
    vertixes_path = Path(vertixes)
    
    mask = vertixes_path.contains_points(all_points)
    ablate_coords=[]
    counter=0
    for i in range(n_pixels):
        for j in range(n_pixels):
            if mask[counter] == False:
                #images_processed[100][i][j]=0
                ablate_coords.append([i,j])
            counter+=1
            
    return ablate_coords

    
def Add_triple_arrays(coords1, coords2, coords3):
    
    coords = []
    if len(coords1)>0: coords.extend(coords1)
    if len(coords2)>0: coords.extend(coords2)
    if len(coords3)>0: coords.extend(coords3)

    #coords = np.vstack(coords)

    #if len(coords1)>0 and len(coords2)>0:
    #    coords = np.array(np.vstack((coords1, coords2)), dtype=np.int16)
    #
    #if len(coords3)>0:
    #    coords = np.array(np.vstack((coords, coords3)), dtype=np.int16)
        

    
    return coords

def Add_double_arrays(bregma_coords_temp, generic_coords):
    coords = np.array(np.vstack((bregma_coords_temp, generic_coords)), dtype=np.int16)
    
    return coords

def Load_areas_and_mask(depth, unit, channel, n_pixels, main_dir, file_dir, file_name, images_aligned, area_names, sides):
    
    #Create set of images, saved coordinates and borders for output
    images_areas = []   #Make list of lists to hold images for each area

    with open(file_dir+file_name+'/depth.txt', "r") as f:
        data = csv.reader(f)
        for row in data:
            depth = row[0]
    
    #Load stimulus evoked areas:
    print "Loading ROIs..."
    for area in area_names:
        for side in sides:
            #print area+" "+side
            area_file = file_dir+ depth+'_' + area+'_'+side
            #print area_file
            if (os.path.exists(area_file+'.npy')==True):
                area_mask = np.load(area_file+'.npy')     #Load coordinates; 2D vectors stacked; Mask out everything else
               
                a = []
                for i in range(len(images_aligned)): #Mask every processed frame in +/- window period
                    a.append(np.ma.array(images_aligned[i], mask=area_mask, fill_value = 0., hard_mask = True))
                                                    
                images_areas.append(a)
    
    return images_areas

def Load_areas_and_mask_old(unit, channel, n_pixels, main_dir, file_dir, file_name, images_processed, area_names, sides):
    
    #Create set of images, saved coordinates and borders for output
    images_areas = []   #Make list of lists to hold images for each area
    coords_save = []    #holds all ablated coords areas - NOT SURE REQUIRED ANYMORE
    borders = []

    with open(file_dir+file_name+'/depth.txt', "r") as f:
        data = csv.reader(f)
        for row in data:
            depth = row[0]
    
    #Load stimulus evoked areas:
    print "Loading ROIs..."
    for area in area_names:
        for side in sides:
            print area+" "+side
            area_file = main_dir + area+'_'+side
            if (os.path.exists(area_file+'.npy')==True):
                area_coords = np.load(area_file+'.npy')     #Load coordinates; 2D vectors stacked; Mask out everything else
                
                #Mask all values = 1.0 (i.e. True) except the area interested in (i.e False)
                tempmask_indexes=np.zeros((n_pixels,n_pixels))+1.0
                for i in range(len(area_coords)):
                    tempmask_indexes[area_coords[i][0]][area_coords[i][1]] = False
                
                a = []
                for i in range(len(images_processed)): #Mask every processed frame in +/- window period
                    a.append(np.ma.array(images_processed[i], mask=tempmask_indexes))
                
                plt.imshow(a[0])
                plt.show()
                #quit()
                
                
                images_areas.append(a)
                
                temp_array=[]
                for i in range(len(tempmask_indexes)):
                    for j in range(len(tempmask_indexes[i])):
                        if tempmask_indexes[i][j]== True: temp_array.append([i,j])

                coords_save.append(temp_array)
    
    print "...done."
    return images_areas, coords_save

def Load_areas_whole_cortex(unit, channel, n_pixels, file_dir, file_name, images_processed):
        
    #Create set of images for each defined area
    images_areas = []   #Make list of lists to hold images for each area
    coords_save = []    #holds all ablated coords areas
    borders = []
    
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
    bregma_coords_temp = bregma_coords.copy()
    for j in range(len(bregma_coords)):  #Remove centreline from all images
        for k in range(n_pixels):
            for l in range(7):
                bregma_coords_temp=np.vstack((bregma_coords_temp, [k,min(n_pixels-1,bregma_coords[1]-3+l)]))
    
    generic_mask_indexes=np.zeros((n_pixels,n_pixels))
    for i in range(len(generic_coords)):
        generic_mask_indexes[generic_coords[i][0]][generic_coords[i][1]] = True
    for i in range(len(bregma_coords_temp)):
        generic_mask_indexes[bregma_coords_temp[i][0]][bregma_coords_temp[i][1]] = True       
    for i in range(len(artifact_coords)):
        generic_mask_indexes[artifact_coords[i][0]][artifact_coords[i][1]] = True       
    
    #Load All Cortex 
    ipsi_coords = []
    ablate_coords = Add_triple_arrays(bregma_coords_temp, generic_coords, artifact_coords)
    
    tempmask_indexes=np.zeros((n_pixels,n_pixels))
    for i in range(len(ablate_coords)):
        tempmask_indexes[ablate_coords[i][0]][ablate_coords[i][1]] = True
    a = []
    for i in range(len(images_processed)):
        a.append(np.ma.array(images_processed[i], mask=tempmask_indexes))
    
    images_areas = a
    coords_save=ablate_coords
    borders= ipsi_coords
    
    return images_areas, coords_save, generic_mask_indexes#, borders

def Compute_static_maps_max(img_rate, window, n_procs, main_dir, file_dir, file_name, n_pixels, unit, channel, n_spikes, ptp):
    ''' Average pre- and post- spike intervals to obtain static map; See also the "max" version of this function'''
    print "Computing static maps"
    
    plot_string = 'all'     #Trigger off all spkes
    #plot_string = '1sec'    #Trigger only on spikes with at least 0.5sec silence on both sides
    #plot_string = 'burst'   #Trigger only off 1st spike in burst, lockout following 0.5 seconds    

    npy_file_name = glob.glob(file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_'+plot_string+'_'+str(window)+'*')
    if len(npy_file_name)==0: return #Skipped cell that doesn't have img_avg; usually the case for cortex recs that have deep cells from subcortex; or zero spike periods
    #print npy_file_name[0]
    
    #npy_file_name = file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_'+plot_string+'_'+str(window)+'sec_window'
    
    images_areas = np.load(npy_file_name[0])

    #Load General mask (removes background)
    generic_mask_file = []
    generic_mask_file = main_dir + 'genericmask.txt'
    if (os.path.exists(generic_mask_file)==True):
        generic_coords = np.int16(np.loadtxt(generic_mask_file))
    
    generic_mask_indexes=np.zeros((n_pixels,n_pixels))
    for i in range(len(generic_coords)):
        generic_mask_indexes[generic_coords[i][0]][generic_coords[i][1]] = True
    #for i in range(len(bregma_coords_temp)):
        #generic_mask_indexes[bregma_coords_temp[i][0]][bregma_coords_temp[i][1]] = True       
    #for i in range(len(artifact_coords)):
        #generic_mask_indexes[artifact_coords[i][0]][artifact_coords[i][1]] = True     

    if True:
        temp_array = []
        for i in range(-int(img_rate), int(img_rate),1):
            temp = np.ma.array(images_areas[int(window*img_rate)+i], mask=generic_mask_indexes, fill_value = 0., hard_mask = True)
            temp_array.append(temp)

        temp_array = np.float32(temp_array)
        images_out2 = np.amax(temp_array,axis=0)
        images_out2 = np.ma.array(images_out2, mask=generic_mask_indexes, fill_value = 0)
        images_out2 = np.ma.filled(images_out2, 0.0)
        images_out2 = np.nan_to_num(images_out2)
        #print "Max: ", np.max(images_out2), np.min(images_out2)
        
        ##images_out3 contains average map from t=0 for folloowing 10 frames (~100-500ms)
        #images_out3 = np.average(temp_array[int(len(temp_array)/2):int(len(temp_array)/2)+10],axis=0)
        #images_out3 = np.ma.array(images_out3, mask=generic_mask_indexes, fill_value = 0)
        #images_out3 = np.ma.filled(images_out3, 0.0)
        #images_out3 = np.nan_to_num(images_out3)
        ##print "Ave: ", np.max(images_out3), np.min(images_out3)
        
        #ax = plt.subplot(1,3,1)
        #plt.imshow(images_out2)
        #ax = plt.subplot(1,3,2)
        #plt.imshow(images_out3)
        #ax = plt.subplot(1,3,3)
        #plt.imshow(images_out3-images_out2)
        #plt.show()

        
        #file_ = file_dir + file_name+'/'+file_name+'_postonly_unit_'+str(unit).zfill(2)+'_channel_'+str(channel).zfill(2)+'_ptp_'+str(ptp).zfill(3)+'_spikes_'+str(n_spikes).zfill(5)
        #images_out2.dump(file_)
        
        np.save(file_dir + file_name+'/'+file_name+'_'+plot_string+'_maxmap_unit_'+str(unit).zfill(2)+'_channel_'+str(channel).zfill(2)+'_ptp_'+str(ptp).zfill(3)+'_spikes_'+str(n_spikes).zfill(5), images_out2)
        #scipy.misc.imsave(file_dir + file_name+'/'+file_name+'_postonly_unit_'+str(unit).zfill(2)+'_channel_'+str(channel).zfill(2)+'_ptp_'+str(ptp).zfill(3)+'_spikes_'+str(n_spikes).zfill(5)+'.png', images_out2)
        #np.save(file_dir + file_name+'/'+file_name+'_'+plot_string+'_avemap_unit_'+str(unit).zfill(2)+'_channel_'+str(channel).zfill(2)+'_ptp_'+str(ptp).zfill(3)+'_spikes_'+str(n_spikes).zfill(5), images_out3)

    #SAVE MIN MAPS ALSO
    if True:
        temp_array = []
        for i in range(-int(img_rate), int(img_rate),1):
            temp = np.ma.array(images_areas[int(window*img_rate)+i], mask=generic_mask_indexes, fill_value = 0., hard_mask = True)
            temp_array.append(temp)

        #plt.imshow(temp_array[0])
        #plt.show()
        temp_array = np.float32(temp_array)
        images_out2 = np.amin(temp_array,axis=0)
        images_out2 = np.ma.array(images_out2, mask=generic_mask_indexes, fill_value = 0.)
        
        images_out2 = np.ma.filled(images_out2, np.min(images_out2))
        
        np.save(file_dir + file_name+'/'+file_name+'_'+plot_string+'_minmap_unit_'+str(unit).zfill(2)+'_channel_'+str(channel).zfill(2)+'_ptp_'+str(ptp).zfill(3)+'_spikes_'+str(n_spikes).zfill(5), images_out2)

    #image_loaded = np.load(file_dir + file_name+'/'+file_name+'_postonly_unit_'+str(unit).zfill(2)+'_channel_'+str(channel).zfill(2)+'_ptp_'+str(ptp).zfill(3)+'_spikes_'+str(n_spikes).zfill(5)+'.npy')

    #print np.max(image_loaded)
    #print np.min(image_loaded)
    #plt.imshow(image_loaded)
    #plt.show()
    #quit()
    
def Load_max_maps(ptps, file_dir, file_name):

    with open(file_dir+file_name+'/depth.txt', "r") as f:
        data = csv.reader(f)
        for row in data:
            depth = row[0]
    with open(file_dir+file_name+'/state.txt', "r") as f:
        data = csv.reader(f)
        for row in data:
            state = row[0]

    max_maps = []
    min_maps = []
    max_maps_nspikes = []
    min_maps_nspikes = []
    max_maps_ptp = []
    min_maps_ptp = []
    max_maps_depth = []
    min_maps_depth = []
    max_maps_state = []
    min_maps_state = []
    max_maps_channel = []
    min_maps_channel = []
    max_maps_unit = []
    counter = 0
    files = os.listdir(file_dir+file_name+'/')
    for file_ in files: #Load all cells from list of experiments
        #print file_
        if ('maxmap' in file_):
        #if ('avemap' in file_):

            n_spikes = int(file_[-9:-4])
            ptp = int(file_[-20:-17])
            channel = int(file_[-27:-25])
            unit = int(file_[-38:-36])
            #if qc[unit]==0: continue
            
            maps = np.load(file_dir+file_name+'/'+file_)
            max_maps.append(maps)
            max_maps_nspikes.append(n_spikes)
            max_maps_ptp.append(ptp)
            max_maps_depth.append(depth)
            max_maps_state.append(state)
            max_maps_channel.append(channel)
            max_maps_unit.append(unit)
            #print channel
            #print "Maxmap spikes: ", n_spikes
            
        if ('minmap' in file_):
            n_spikes = int(file_[-9:-4])
            ptp = int(file_[-20:-17])
            channel = int(file_[-27:-25])

            maps = np.load(file_dir+file_name+'/'+file_)
            min_maps.append(maps)
            min_maps_nspikes.append(n_spikes)
            min_maps_ptp.append(ptp)
            min_maps_depth.append(depth)
            min_maps_state.append(state)
            min_maps_channel.append(channel)
            #print "Minmap spikes: ", n_spikes
    #max_maps = max_maps/counter/255
    
    return max_maps, min_maps, max_maps_nspikes, min_maps_nspikes, max_maps_ptp, min_maps_ptp, max_maps_depth,min_maps_depth,max_maps_state,min_maps_state,max_maps_channel,min_maps_channel,max_maps_unit


def Mp_max((temp_array1)):
    global coords_temp
    
    for j in range(len(coords_temp)):
        temp_array1[coords_temp[j][0]][coords_temp[j][1]] = -1000

    return temp_array1

    
def Mp_min((temp_array2)):
    global coords_temp

    for j in range(len(coords_temp)):
        temp_array2[coords_temp[j][0]][coords_temp[j][1]] = 1000

    return temp_array2

def Zero_images((images_areas)):
    global coords_temp
    
    for i in range(len(images_areas)):
        print "Processing time: ", i
        for j in range(len(coords_temp)):
            print "Processing coordinate: ", j
            images_areas[i][min(255,int(coords_temp[j][0]))][min(255,int(coords_temp[j][1]))]=0
            
    return images_areas

def Search_max_min(images_processed, file_dir, file_name, img_rate, window, n_procs, area_names, depth, sides, generic_mask_indexes):

    global coords_temp

    Max_plot = []
    Min_plot = []
    Max_index = []
    Min_index = []
    Max_pixel_value = []
    Min_pixel_value = []
    Images_areas = []

    #Search max/min for each ROI defined area
    counter=0
    for area in area_names:
        for side in sides:
            print "Searching max/min for area # ", area, " ", side
            
            area_file = file_dir+ depth+'_' + area+'_'+side

            if (os.path.exists(area_file+'.npy')==True):
                area_mask = np.load(area_file+'.npy')     #Load coordinates; 2D vectors stacked; Mask out everything else
            
            coords_temp = area_mask

            temp_array = np.ma.array(np.zeros(images_processed.shape), mask=True)
            for i in range(len(images_processed)):
                temp_array[i] = np.ma.array(images_processed[i], mask=coords_temp, fill_value = 0., hard_mask = True)
            
            #ax =plt.subplot(131)
            #v_min = np.nanmin(images_processed[100])
            #v_max = np.nanmax(images_processed[100])
            #plt.imshow(images_processed[100], vmin=v_min, vmax=v_max)

            #ax =plt.subplot(132)
            #plt.imshow(coords_temp)
            
            #ax =plt.subplot(133)
            #plt.imshow(temp_array[100])
            #plt.show()
            
            temp_max_array = []
            temp_min_array = []
            temp_array1 = np.ma.array(np.zeros((int(img_rate*2),256,256)), mask=True)
            temp_array2 = np.ma.array(np.zeros((int(img_rate*2),256,256)), mask=True)
            for i in range(int(img_rate*2)):     #This searches +/- 1 sec from time = 0 sec; OR Maximum to begining of data;
                temp_array1[i] = np.ma.array(temp_array[int(window*img_rate)-int(img_rate)+i]).copy()
                temp_array2[i] = np.ma.array(temp_array[int(window*img_rate)-int(img_rate)+i]).copy()

            temp_array1._sharedmask=False
            temp_array2._sharedmask=False


            #Set masked background values to very large/very low values to not confuse search for max/min values
            pool = mp.Pool(n_procs)
            temp_max_array.extend(pool.map(Mp_max, temp_array1))
            pool.close()

            pool = mp.Pool(n_procs)
            temp_min_array.extend(pool.map(Mp_min, temp_array2))
            pool.close()

            #Search for max value using 1D unravel of temp_arrays; assign location of max/min index; detect max/min values overall;
            temp1_array = np.ma.array(np.zeros((int(img_rate*2),256,256)), mask=True)
            for k in range(len(temp_max_array)):
                temp1_array[k] = np.ma.array(temp_max_array[k], mask=area_mask) #np.ma.array(temp_max_array)

            temp_max_array = temp_array1
            max_index = np.unravel_index(np.argmax(temp_max_array), temp_max_array.shape)
            max_pixel_value = temp_array[int(window*img_rate)-int(img_rate)+max_index[0]][max_index[1]][max_index[2]]

            #temp1_array = np.ma.array(np.zeros((img_rate*2,256,256)), mask=True)
            #for k in range(len(temp_max_array)):
                #temp1_array[k] = np.ma.array(temp_min_array[k], mask=area_mask) #np.ma.array(temp_max_array)
            temp_min_array = temp_array2
            min_index = np.unravel_index(np.argmin(temp_min_array), temp_min_array.shape)
            min_pixel_value = temp_array[int(window*img_rate)-int(img_rate)+min_index[0]][min_index[1]][min_index[2]]
            
            #max/min_index at index [0] contain time slice info; needed for max/min_pixel_value above, but not after;
            max_index = max_index[1:] #Reduce back to 2D arrays
            min_index = min_index[1:]
           
            max_plot = []
            min_plot = []
            
            for i in range(len(temp_array)):
                max_plot.append(temp_array[i][max_index[0]][max_index[1]])
                min_plot.append(temp_array[i][min_index[0]][min_index[1]])
              
            Max_plot.append(max_plot)
            Min_plot.append(min_plot)
            Max_index.append(max_index)
            Min_index.append(min_index)
            Max_pixel_value.append(max_pixel_value)
            Min_pixel_value.append(min_pixel_value)
            Images_areas.append(temp_array)
            
            counter+=1


    return Max_plot, Min_plot, Max_pixel_value, Min_pixel_value, Max_index, Min_index, area_names, Images_areas

def Average_roi(images_areas, img_rate, window, n_procs, area_names, sides):
    
    global coords_temp

    ave_plot = []
    print "Averaging areas..."
    
    #Search max/min for each ROI defined area
    counter=0
    for area in area_names:
        for side in sides:
            
            max_plot = []
            for i in range(len(images_areas[counter])):     #This searches +/- 1 sec from time = 0 sec; OR Maximum to begining of data;
                value = np.nanmean(images_areas[counter][i])
                max_plot.append(value)
            
            ave_plot.append(max_plot)

            counter+=1

    return ave_plot
    


def Search_max_min_single_area(images_areas, img_rate, window, coords_save, n_procs, tau):

    global coords_temp

    #area_names = [ 'hindlimb', 'barrel', 'motor', 'visual', 'retrosplenial', 'acc', 'allcortex']

    Max_plot = []
    Min_plot = []
    Max_index = []
    Min_index = []
    Max_pixel_value = []
    Min_pixel_value = []

    #Loop over each cortical area
    print "Searching max/min for single area map"
    coords_temp = coords_save
    #print coords_temp

    #print len(tau)
    #print len(images_areas)

    temp_max_array = []
    temp_min_array = []
    temp_array1 = []
    temp_array2 = []

    for i in range(60):     #This searches +/- 1 sec from time = 0 sec; OR Maximum to begning of data;
        temp_array1.append(np.array(images_areas[60+i]).copy())
        temp_array2.append(np.array(images_areas[60+i]).copy())

    #Set masked background values to very large/very low values to not confuse search for max/min values
    pool = mp.Pool(n_procs)
    temp_max_array.extend(pool.map(Mp_max, temp_array1))
    pool.close()

    pool = mp.Pool(n_procs)
    temp_min_array.extend(pool.map(Mp_min, temp_array2))
    pool.close()

    temp_max_array = np.array(temp_max_array)
    max_index = np.unravel_index(np.argmax(temp_max_array), temp_max_array.shape)
    max_pixel_value = images_areas[60+max_index[0]][max_index[1]][max_index[2]]

    temp_min_array = np.array(temp_min_array)
    min_index = np.unravel_index(np.argmin(temp_min_array), temp_min_array.shape)
    min_pixel_value = images_areas[60+min_index[0]][min_index[1]][min_index[2]]

    #print max_index
    #plt.imshow(images_areas[60+max_index[0]], vmin=min_pixel_value,vmax=max_pixel_value)
    #plt.show()
    #quit()
    
    max_index = max_index[1:] #Reduce back to 2D arrays
    min_index = min_index[1:]
   
    max_plot = []
    min_plot = []
    
    for i in range(len(images_areas)):
        #Time course curve data - before zeroing out the pixels
        max_plot.append(images_areas[i][max_index[0]][max_index[1]])
        min_plot.append(images_areas[i][min_index[0]][min_index[1]])
      
    Max_plot.append(max_plot)
    Min_plot.append(min_plot)
    Max_index.append(max_index)
    Min_index.append(min_index)
    Max_pixel_value.append(max_pixel_value)
    Min_pixel_value.append(min_pixel_value)
    
    return Max_plot, Min_plot, Max_pixel_value, Min_pixel_value, Max_index, Min_index


def Save_time_course(unit, channel, spikes, Max_plot, Min_plot, Max_index, Min_index, window, len_frame, file_dir, file_name, area_names, sides, plot_string):
    
    #colors = ['blue', 'green', 'orange', 'brown', 'red', 'magenta', 'black', 'yellow']

    #temp_array saves data to file; first line contains number of spikes and area names (though have to be careful with area_names as they  
    temp_array = []
    temp_array.append([len(spikes)])  #Save # spikes within imaging period
    temp_array2 = []
    for area in area_names:
        for side in sides:
            temp_array2.append(area+"_"+side)   #Save names of areas recorded
    #temp_array2.append('allcortex')
    #temp_array.append(temp_array2)
    
    for i in range(len(Max_plot)):
        ax = plt.subplot(8,2,i+1)
        max_plot=np.array(Max_plot[i])
        min_plot=np.array(Min_plot[i])

        temp_array.append(min_plot)
        temp_array.append(max_plot)

        xx = np.arange(-window, window, len_frame)
        xx = xx[0: len(max_plot)]

        ax.plot(xx, max_plot, color='black', linewidth=2)
        ax.plot(xx, min_plot, color='blue', linewidth=2)
        ax.plot([-3,3],[0,0], color='black')
        ax.plot([0,0],[min(min_plot),max(max_plot)], color='black')
        ax.set_ylim(-.04, .041)
        ax.tick_params(axis='both', which='major', labelsize=8)
        ax.yaxis.set_ticks(np.arange(-0.04,0.041,0.04))
        
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        if i==0: 
            plt.title("Left")
            ax.get_yaxis().set_visible(True)
            plt.ylabel(area_names[int(i/2)], fontsize=8)
        elif (i%2)==0 and i<14:
            ax.get_yaxis().set_visible(True)
            plt.ylabel(area_names[int(i/2)], fontsize=8)
        elif i==1:
            #ax.get_xaxis().set_visible(True)
            plt.title("Right")
        elif i==13:
            ax.get_xaxis().set_visible(True)
        #elif i==14:
        #   ax.get_yaxis().set_visible(True)
        #    plt.ylabel("allcortex", fontsize=8)

    plt.suptitle(file_name + " Unit: " +str(unit).zfill(2) + " Channel: " + str(channel).zfill(2)+
    " No. of spikes: "+ str(len(spikes)))

    #npy_file_name = file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_'+plot_string+'_'+str(window)+'sec_window'

    plt.savefig(file_dir + file_name + '/time_course_plot_' + file_name+'_'+plot_string+ '_'+str(window)+'sec_window_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'.png', fontsize = 20)
    plt.close()

    #Add max_index and min_index information to the end of the time_course_data*.txt file
    for i in range(len(Max_index)):
        temp_array.append(Max_index[i])
        temp_array.append(Min_index[i])

    #Save time_course_data
    with open(file_dir + file_name + '/time_course_data_' + file_name+'_'+plot_string+ '_'+str(window)+'sec_window_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+"_"+area_names[0], "w") as f:
        writer = csv.writer(f)
        writer.writerows(temp_array)


def Plot_time_course(unit, channel, spikes, Max_plot, Min_plot, window, len_frame, file_dir, file_name, area_names):
    print "Plotting time course from file"
    
    with open(file_dir + file_name + '/time_course_data_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2), "r") as f:
        data = csv.reader(f)
        for row in data:
            print row

def Plot_matrix_maps(average_areas, file_dir, file_name, area_names, img_rate, unit, spikes, channel, ptp):
    
    #fig = plt.figure()
    plt.close()
    gs = gridspec.GridSpec(2,6)
    
    fig, ax = plt.subplots(nrows=1, ncols=2)

    #Compute time courses
    
    img1 = []
    for k in range(len(area_names)):
       img1.append(np.float16(average_areas[k*2])*1E2)
    v_abs1 = max(np.max(np.array(img1)), -np.min(np.array(img1)))

    img2 = []
    for k in range(len(area_names)):
       img2.append(np.float32(average_areas[k*2+1])*1E2)
    v_abs2 = max(np.max(np.array(img2)), - np.min(np.array(img2)))

    global_max = max(v_abs1,v_abs2)

    #****** Plot left hemisphere
    #ax = plt.subplot(gs[0,0:1])
    ax1 = plt.subplot(121)

    im  = ax1.imshow(img1, aspect='auto', cmap=plt.get_cmap('jet'), vmin=-global_max, vmax=global_max, interpolation='none')
    
    yy = np.arange(0,len(area_names),1)
    #labels = [ 'hindlimb', 'forelimb', 'barrel', 'motor', 'visual', 'retrosplenial', 'acc']
    plt.yticks(yy, area_names) 
    plt.tick_params(axis='both', which='both', labelsize=8)

    original_xx = np.arange(0,len(average_areas[0])+2,img_rate)
    xx = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
    plt.ylim(-0.5,-0.5+len(area_names))
    plt.xlim(0,len(average_areas[0]))
    plt.xticks(original_xx, xx)
    
    plt.plot([3*img_rate,3*img_rate],[-0.5,-0.5+len(area_names)], color='black', linewidth=2)
    for i in range(len(area_names)):
        plt.plot([0,len(average_areas[0])+2],[-0.5+i,-0.5+i], 'r--', color='black', linewidth=1)
    for i in range(7):
        plt.plot([0+i*img_rate,0+i*img_rate],[-0.5,-0.5+len(area_names)], 'r--', color='black', linewidth=1)
    plt.title("Left (DF/F max: "+str(round(np.max(img1),2))+"  min: "+str(round(np.min(img1),2))+")", fontsize=10) 

    np.save(file_dir+file_name+'/'+file_name+'_matrix_unit_'+str(unit).zfill(2)+'_left', img1)
    
    
    #****** Plot right hemisphere
    #ax = plt.subplot(gs[0,1:2])
    ax2 = plt.subplot(122)
    im = ax2.imshow(img2, aspect='auto', cmap=plt.get_cmap('jet'), vmin=-global_max, vmax=global_max, interpolation='none')
    plt.tick_params(axis='both', which='both', labelsize=8)
    ax2.get_yaxis().set_visible(False)

    plt.ylim(-0.5,-0.5+len(area_names))
    plt.xlim(0,len(average_areas[0]))
    plt.xticks(original_xx, xx)
    plt.plot([3*img_rate,3*img_rate],[-0.5,-0.5+len(area_names)], color='black', linewidth=2)
    for i in range(len(area_names)):
        plt.plot([0,len(average_areas[0])+2],[-0.5+i,-0.5+i], 'r--', color='black', linewidth=1)
    for i in range(7):
        plt.plot([0+i*img_rate,0+i*img_rate],[-0.5,-0.5+len(area_names)], 'r--', color='black', linewidth=1)

    plt.title("Right (DF/F max: "+str(round(np.max(img2),2))+"  min: "+str(round(np.min(img2),2))+")", fontsize=10) 


    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax)
    plt.suptitle(file_name + "  unit: " + str(unit)+ "  ptp: " + str(ptp)+ "uV  spikes: "+str(len(spikes))+ "  ch: "+str(channel))

    np.save(file_dir+file_name+'/'+file_name+'_matrix_unit_'+str(unit).zfill(2)+'_right', img2)


    if (ptp>35) and (len(spikes)>128) and (global_max>0.5):
        plt.savefig(file_dir + file_name + '/'+file_name+'_matrix_unit_'+str(unit).zfill(2)+'.png', dpi=100)#, fontsize = 20)
        plt.close()
    else:
        plt.close()
        return
        
    
    return
    
    quit()

    #****** Plot STMTD - traces
    ax = plt.subplot(gs[0,2:4])
    print np.array(average_areas).shape
    ax.get_yaxis().set_visible(False)
    
    original_xx = np.arange(0,len(average_areas[0])+2,img_rate)
    xx = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
    plt.xlim(0,len(average_areas[0]))
    plt.ylim(-0.015,0.13)
    plt.xticks(original_xx, xx)
    plt.plot([img_rate*3,img_rate*3],[-0.02,0.15], 'r--', color='black', linewidth=2)
    
    plt.title("STMTD")
    for k in range(7):
        plt.plot(np.array(average_areas[k*2])+k*0.02)
        plt.plot([0,img_rate*6],[k*0.02,k*0.02], 'r--', color='black')
       
    #****** Plot Maxmaps
    ax = plt.subplot(gs[0,4:6])
    map_name = glob.glob(file_dir + file_name+'/'+file_name+'_maxmap_unit_'+str(unit).zfill(2)+'*.npy')
    #print "MAP: ", map_name[0]
    #quit()
    img=np.load(map_name[0])
    
    v_max = 0
    v_min = np.min(img)
    for k in range(len(img)):
        temp = np.max(img[k])
        if v_max<temp: v_max = temp
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_visible(False)

    plt.title("STM  -  Maxpixel Method (-1..+1 sec)\n DF/F min: "+str(round(v_min*100.,2)) + "  max: "+str(round(v_max*100.,2)))
    plt.imshow(img, vmin=v_min, vmax=v_max)
                        
    #******* Plot static map
    map_name = glob.glob(file_dir + file_name+'/img_avg_'+file_name+'_unit'+str(unit).zfill(2)+'*.npy')
    img_data=np.load(map_name[0])
    v_max = np.max(img_data)
    v_min = np.min(img_data)

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
                #print [k,min(n_pixels-1,bregma_coords[j][1]-3+l)]
                bregma_coords_temp.append([k,min(n_pixels-1,bregma_coords[j][1]-3+l)])
    
    generic_mask_indexes=np.zeros((n_pixels,n_pixels))
    for i in range(len(generic_coords)):
        generic_mask_indexes[generic_coords[i][0]][generic_coords[i][1]] = True
    for i in range(len(bregma_coords_temp)):
        generic_mask_indexes[bregma_coords_temp[i][0]][bregma_coords_temp[i][1]] = True       
    for i in range(len(artifact_coords)):
        generic_mask_indexes[artifact_coords[i][0]][artifact_coords[i][1]] = True     

    img_plot = []
    for t in range(0,6,1):
        print "time index: ", t
        temp_array = []
        for i in range(t*int(img_rate), (t+1)*int(img_rate),1):
            temp = np.ma.array(img_data[i], mask=generic_mask_indexes, fill_value = 0., hard_mask = True)
            temp_array.append(temp)

        temp_array = np.float32(temp_array)
        
        if False:
            #max OR min map
            images_out2 = np.amax(temp_array,axis=0)
            #images_out2 = np.amin(temp_array,axis=0)
            images_out2 = np.ma.array(images_out2, mask=generic_mask_indexes, fill_value = 0)
            images_out2 = np.ma.filled(images_out2, 0.0)
            images_out2 = np.nan_to_num(images_out2)
        
        if True:
            #average map 
            images_out2 = np.average(temp_array,axis=0)
            images_out2 = np.ma.array(images_out2, mask=generic_mask_indexes, fill_value = 0)
            images_out2 = np.ma.filled(images_out2, 0.0)
            images_out2 = np.nan_to_num(images_out2)
        
        img_plot.append(images_out2)
    
    times = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
    v_min=1.0
    v_max=-1.0
    for p in range(6):
        v_min = min(v_min,np.min(img_plot[p]))
        v_max = max(v_max,np.max(img_plot[p]))
        
    for k in range(6):
        ax = plt.subplot(gs[1,k])
        plt.title(str(times[k])+" .. "+str(times[k+1])+" sec")
        #plt.imshow(img_plot[k], vmin=v_min, vmax=v_max)
        plt.imshow(img_plot[k], vmin=np.min(img_plot[k]), vmax=np.max(img_plot[k]))
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)


    plt.suptitle(file_name+"\n  unit: "+str(unit)+ "  n_spikes: "+str(n_spikes), fontsize=20)

    plt.show()
    quit()

def plot_ROI(roi_left,roi_right, unit, area_names, file_name):


    fig, ax = plt.subplots(nrows=1, ncols=2)

    #Plot left hemisphere
    ax1 = plt.subplot(121)
    plt.imshow(roi_left, aspect='auto', vmin = -global_max, vmax=global_max, cmap=plt.get_cmap('jet'), interpolation='none')
    plt.title("DF/F min: "+ str(round(np.min(roi_left),2))+ "  max: "+str(round(np.max(roi_left),2)))

    img_rate = len(roi_left[0])/6
    original_xx = np.arange(0,len(roi_left[0])+2,img_rate)
    xx = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
    plt.xticks(original_xx, xx)
    plt.xlim(0,len(roi_left[0]))

    yy = np.arange(0,len(area_names),1)
    plt.yticks(yy, area_names) 
    plt.ylim(-0.5,-0.5+len(area_names))

    plt.tick_params(axis='both', which='both', labelsize=15)

    plt.plot([3*img_rate,3*img_rate],[-0.5,-0.5+len(area_names)], color='black', linewidth=2)
    for i in range(len(area_names)):
        plt.plot([0,len(roi_left[0])+2],[-0.5+i,-0.5+i], 'r--', color='black', linewidth=1)
    for i in range(7):
        plt.plot([0+i*img_rate,0+i*img_rate],[-0.5,-0.5+len(area_names)], 'r--', color='black', linewidth=1)

    #Plot right hemisphere
    ax2 = plt.subplot(122)
    im = ax2.imshow(roi_right, aspect='auto', vmin = -global_max, vmax=global_max, cmap=plt.get_cmap('jet'), interpolation='none')
    plt.title("DF/F min: "+ str(round(np.min(roi_right),2))+ "  max: "+str(round(np.max(roi_right),2)))

    #Plot 
    img_rate = len(roi_left[0])/6
    original_xx = np.arange(0,len(roi_left[0])+2,img_rate)
    xx = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
    plt.xticks(original_xx, xx)
    plt.xlim(0,len(roi_left[0]))
    
    ax2.get_yaxis().set_visible(False)
    plt.ylim(-0.5,-0.5+len(area_names))

    plt.tick_params(axis='both', which='both', labelsize=15)

    plt.plot([3*img_rate,3*img_rate],[-0.5,-0.5+len(area_names)], color='black', linewidth=2)
    for i in range(len(area_names)):
        plt.plot([0,len(roi_left[0])+2],[-0.5+i,-0.5+i], 'r--', color='black', linewidth=1)
    for i in range(7):
        plt.plot([0+i*img_rate,0+i*img_rate],[-0.5,-0.5+len(area_names)], 'r--', color='black', linewidth=1)

    #Plot color bar
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.1, 0.05, 0.8])
    fig.colorbar(im, cax=cbar_ax)

    #figManager = plt.get_current_fig_manager()
    #figManager.window.showMaximized()
    plt.suptitle(file_name + "  unit: " + str(unit), fontsize=20)
    plt.show()

def Animate_images(unit, channel, window, img_rate, Images_areas, file_dir, file_name, n_pixels, spikes, plot_string, n_procs, generic_mask_indexes, 
    Max_plot, Min_plot, Max_pixel_value, Min_pixel_value, Max_index, Min_index, area_names, sides, depth):
    
    print "Generating Animation. No. of Frames: ", len(Images_areas[0])
    
    colors=['blue','red', 'green','violet','lightseagreen','lightsalmon','indianred','lightsalmon','pink','darkolivegreen','blue','red', 'green','violet','lightseagreen','lightsalmon','indianred','lightsalmon','pink','darkolivegreen']

    #***************** Process plots *****************
    plot_frames = []
    x = np.arange(-window, window, 1./img_rate)
    x = x[:len(Max_plot[0])]
    Max_plot = np.array(Max_plot)

    print "Number of area recordings: ", len(Max_plot)
    ranges = np.arange(0,len(x),1)

    #Divide data into nchunks for some reason works faster; TO DO: PARALLELIZE
    n_chunks = 8
    range_chunks = np.array_split(ranges,n_chunks)
    
    plt.close()
    
    #*********** PLOT CURVE GRAPHS ************************
    if True: #PLOT ALL CURVES FOR ALL AREAS
        print "... processing plot frames..."
        for chunk in range(n_chunks):
            fig = plt.figure()
            fig.add_subplot(111)
            fig.set_size_inches(5, 20)
            #fig.canvas.draw()
            
            #FIXED LINES FOR EACH PLOT
            for i in range(len(Max_plot)+2): #need to make figs at top of graph also                
                plt.plot(x, [int(i/2)*0.10]*len(x), color='black', linewidth = 2)
                plt.plot(x, [int(i/2)*0.10+0.03]*len(x), 'r--', color=colors[int(i/2)], alpha=0.5)
                plt.plot(x, [int(i/2)*0.10-0.03]*len(x), 'r--', color=colors[int(i/2)], alpha=0.5)
                plt.axhspan(int(i/2)*0.10+0.03, int(i/2)*0.10-0.03, color='black', alpha=0.2, lw=0)   
            plt.xlabel("Time (s)")
            plt.ylim(-0.06, 0.90)
            plt.plot([-window,window],[-100,100], color='black')

            plt.tick_params(
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            labelbottom='off') # 

            #DRAW CURVES OVER TIME
            for k in range_chunks[chunk]:
                for i in range(len(Max_plot)):

                    plt.plot(x[:k], Max_plot[i][:k]+(int(i/2)+1)*0.10, color=colors[i], linewidth=2) #Plot individual overlayed left/rigth curves
                    
                fig.savefig(file_dir + file_name+'/figs_'+str(k).zfill(3)+'.jpg', dpi=40)
                #fig.canvas.draw()
                #fig.canvas.flush_events()
                #data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
                #data = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))
                #scipy.misc.toimage(data).save(file_dir + file_name+'/figs_'+str(k).zfill(3)+'.jpg')
            plt.close(fig)

        #Make vid1 containing time course curves
        devnull = open(os.devnull, 'w')
        subprocess.call("ffmpeg -f image2 -r 15 -i " + file_dir+file_name+ "/figs_%03d.jpg  -vcodec libx264 -y -vf scale=480:960 "+file_dir+file_name+'/vid1.mp4', shell=True,stdout=devnull, stderr=devnull)
        os.system('rm '+file_dir+file_name+'/figs*')


    #***************** Process images ****************
    Images_areas = np.array(Images_areas)
    Images_areas = np.swapaxes(Images_areas, 0, 1)
    
    #Generate borders around recorded areas and mask them out in order for color to work properly
    img=Image.new("RGBA", (n_pixels,n_pixels),(0,0,0))
    if True:
        if True:
            draw = ImageDraw.Draw(img)
            counter=0
            for area in area_names:
                for side in sides:
                    area_file = file_dir+ depth+'_' + area+'_'+side+'_contour'
                    if (os.path.exists(area_file+'.npy')==True):
                        area_contour = np.load(area_file+'.npy')     #Load coordinates; 2D vectors stacked; Mask out everything else

                    for i in range(len(area_contour)-1):  #Skip last border which represents all cortex
                        draw.line((area_contour[i][1],area_contour[i][0],area_contour[i+1][1],area_contour[i+1][0]), fill = colors[counter], width = 1)
                counter+=1
        borders_=asarray(img)

        if True:
            #Generate borders mask to zero out [Ca] image pixels to get color right
            borders_mask = np.zeros((n_pixels,n_pixels),dtype=bool)
            for i in range(n_pixels):
                for j in range(n_pixels):
                    if np.sum(borders_[i][j])>255:
                        borders_mask[i][j]=True

        #if False:
            #draw = ImageDraw.Draw(img)
            #for area in area_names:
                #for side in sides:
                    #area_file = file_dir+ depth+'_' + area+'_'+side
                    #if (os.path.exists(area_file+'.npy')==True):
                        #area_mask = np.load(area_file+'.npy')     #Load coordinates; 2D vectors stacked; Mask out everything else

                    ##for i in range(len(area_mask)):  #Skip last border which represents all cortex
                    ##    for j in range(len(area_mask[[0]])):
                    ##        if area_mask[i][j]==True:
                    ##            print i,j
                    ##            draw.point((i,j), fill = 'white')

                
    #Generate time index text and convert it to array for adding to frames
    font = ImageFont.truetype("/usr/share/fonts/truetype/ttf-dejavu/DejaVuSansMono.ttf",9)
    time_text=[]
    for i in range(len(Images_areas)):
        img=Image.new("RGBA", (n_pixels,n_pixels),(0,0,0))
        draw = ImageDraw.Draw(img)
        #if n_pixels>128:
        #    draw.text((0, n_pixels-10),"Spikes: " + str(len(spikes))+ " Time: "+str(format((float(i)/len(Images_areas)-.5)*6.0,'.2f')),(255,255,255),font=font)
        #else:
        draw.text((10, 10),"Spikes: " + str(len(spikes))+ " Time: "+str(format((float(i)/len(Images_areas)-.5)*6.0,'.2f')),(255,255,255),font=font)
        draw = ImageDraw.Draw(img)
        time_text.append(asarray(img))

    label_text=[]
    for i in range(len(Images_areas)):
        img=Image.new("RGBA", (n_pixels,n_pixels),(0,0,0))
        draw = ImageDraw.Draw(img)
        #if n_pixels>128:
        #    draw.text((0, n_pixels-20),file_name + " Unit: "+str(unit),(255,255,255),font=font)
        #else:
        draw.text((10, 0),file_name + " Unit: "+str(unit),(255,255,255),font=font)
        draw = ImageDraw.Draw(img)
        label_text.append(asarray(img))
    
    #Generate image frames as arrays and save to .pngs
    images_out = []
    my_cmap = matplotlib.cm.get_cmap('jet')
    v_min = min(Min_pixel_value)
    v_max = max(Max_pixel_value)

    print "... processing images frames..."
    for i in range(len(Images_areas)):
        temp_img = np.float32(Images_areas[i][0])
        temp_img = ndimage.gaussian_filter(temp_img, sigma=.5)
        temp_img = (temp_img - v_min)/(v_max-v_min)
        masked_data = np.ma.array(temp_img, mask=generic_mask_indexes)
        
        if True: masked_data = np.ma.array(masked_data, mask=borders_mask)
        
        images_out = my_cmap(masked_data, bytes=True) #, vmin=min(Min_pixel_v), vmax=max(Max_pixel_v), origin='lower') #cmap=my_cmap, clim=[0.9, 1]) #, cmap=cmap, interpolation='nearest')
        
        if True: images_out = borders_+ images_out + time_text[i] + label_text[i]

        scipy.misc.imsave(file_dir + file_name+'/figs_'+str(i).zfill(3)+'.png', images_out)
        #scipy.misc.toimage(images_out, cmin=v_min, cmax=v_max).save(file_dir + file_name+'/figs_'+str(i).zfill(3)+'.png')
   
    print "... making vids ..."

    #Make video 2 - [Ca] imaging vid
    devnull = open(os.devnull, 'w')
    subprocess.call("ffmpeg -f image2 -r 15 -i " + file_dir+file_name+ "/figs_%03d.png  -vcodec libx264 -y -vf scale=960:960 "+file_dir+file_name+'/vid2.mp4', shell=True, stdout=devnull, stderr=devnull)
    os.system('rm '+file_dir+file_name+'/figs*')

    #Combine videos into 960 x 480 vid
    devnull = open(os.devnull, 'w')
    subprocess.call("ffmpeg -i "+ file_dir+file_name+"/vid2.mp4 -i " + file_dir+file_name+"/vid1.mp4 -filter_complex '[0:v]pad=iw*2:ih[int];[int][1:v]overlay=W/2:0[vid]' -map [vid] -c:v libx264 -crf 23 -preset veryfast -y "
    + file_dir+file_name+'/'+file_name+'_unit'+str(unit)+'_'+plot_string+'_'+str(window)+'sec_window1.mp4', shell=True,stdout=devnull, stderr=devnull)

    #Crop video
    devnull = open(os.devnull, 'w')
    subprocess.call("ffmpeg -i " + file_dir+file_name+'/'+file_name+'_unit'+str(unit)+'_'+plot_string+'_'+str(window)+"sec_window1.mp4 -filter:v 'crop=1440:960:0:0' -c:a copy -y " 
    + file_dir+file_name+'/'+file_name+'_unit'+str(unit)+'_'+plot_string+'_'+str(window)+'sec_window.mp4', shell=True,stdout=devnull, stderr=devnull)

    #Delete old videos
    subprocess.call("rm " + file_dir + file_name+ "/vid1.mp4 "+ file_dir + file_name+ "/vid2.mp4 "+ file_dir+file_name+'/'+file_name+'_unit'+str(unit)+'_'+plot_string+'_'+str(window)+'sec_window1.mp4', shell=True, stdout=devnull, stderr=devnull)

    print "Finished unit: ", unit
    print ""
    print ""
    print ""
    #quit()

def Animate_images_decorrelated(unit, channel, window, Images_areas, file_dir, file_name, n_pixels, spikes, plot_string, n_procs, generic_mask_indexes,Max_pixel_value, Min_pixel_value):

    print "Generating Animation. No. of Frames: ", len(Images_areas)
    
    #Load Bregma coords for computing contralateral areas
    bregma_mask_file = file_dir + 'bregmamask.txt'
    bregma_coords = np.loadtxt(bregma_mask_file)
    
    colors = ['blue', 'green', 'orange', 'brown', 'red', 'magenta', 'black']

    #Prepare images;
    Images_areas = np.array(Images_areas)
    #print len(Images_areas)
    
    #plt.imshow(Images_areas[0])
    #plt.show()
    #quit()
    #Images_areas = np.swapaxes(Images_areas, 0, 1)

    #Generate time index text and convert it to array for adding to frames
    font = ImageFont.truetype("/usr/share/fonts/truetype/ttf-dejavu/DejaVuSansMono.ttf",9)
    time_text=[]
    for i in range(len(Images_areas)):
        img=Image.new("RGBA", (n_pixels,n_pixels),(0,0,0))
        draw = ImageDraw.Draw(img)
        if n_pixels>128:
            draw.text((0, n_pixels-10),"Spikes: " + str(len(spikes))+ " Time: "+str(format((float(i)/len(Images_areas)-.5)*6.0,'.2f')),(255,255,255),font=font)
        else:
            draw.text((10, 10),"Spikes: " + str(len(spikes))+ " Time: "+str(format((float(i)/len(Images_areas)-.5)*6.0,'.2f')),(255,255,255),font=font)
        draw = ImageDraw.Draw(img)
        time_text.append(asarray(img))

    label_text=[]
    for i in range(len(Images_areas)):
        img=Image.new("RGBA", (n_pixels,n_pixels),(0,0,0))
        draw = ImageDraw.Draw(img)
        if n_pixels>128:
            draw.text((0, n_pixels-20),file_name + " Unit: "+str(unit),(255,255,255),font=font)
        else:
            draw.text((10, 0),file_name + " Unit: "+str(unit),(255,255,255),font=font)
        draw = ImageDraw.Draw(img)
        label_text.append(asarray(img))

    #Generate image frames as arrays and save to .pngs
    images_out = []
    my_cmap = matplotlib.cm.get_cmap('jet')
    v_min = min(Min_pixel_value)
    v_max = max(Max_pixel_value)
    print "v_min: ", v_min
    print "v_max: ", v_max
    #v_min = 0
    #v_max = .5
    
    for i in range(len(Images_areas)):
        #print "Frame: ", i
        #temp_img = np.float32(Images_areas[i])
        #temp_img = ndimage.gaussian_filter(temp_img, sigma=.5)
        temp_img = (Images_areas[i] - v_min)/(v_max-v_min)
        masked_data = np.ma.array(temp_img, mask=generic_mask_indexes)
        images_out = my_cmap(masked_data, bytes=True) #, vmin=min(Min_pixel_v), vmax=max(Max_pixel_v), origin='lower') #cmap=my_cmap, clim=[0.9, 1]) #, cmap=cmap, interpolation='nearest')
        images_out = images_out + time_text[i] + label_text[i]
        scipy.misc.imsave(file_dir + file_name+'/figs_'+str(i).zfill(3)+'.png', images_out)#, vmin=v_min, vmax=v_max)
        #scipy.misc.toimage(images_out, low=v_min, high=v_max).save(file_dir + file_name+'/figs_'+str(i).zfill(3)+'.png')

    #Combine .png files and exit
    import subprocess
    devnull = open(os.devnull, 'w')
    subprocess.call("ffmpeg -f image2 -r 15 -i " + file_dir+file_name+ "/figs_%03d.png  -vcodec libx264 -y "+file_dir+file_name+'/'+file_name+
    '_unit'+str(unit)+'_'+plot_string+'_'+str(window)+'sec_window_decorrelated.mp4', shell=True,stdout=devnull, stderr=devnull)
    os.system('rm '+file_dir+file_name+'/figs*')
        
    print "Finished unit: ", unit
    print ""
    print ""
    quit()


def lfp_img_correlation(file_dir, file_name, images_raw,img_start,img_end,len_frame,img_rate,n_pixels,n_frames, img_times):

    #Load Generic Mask
    coords_generic=[]
    if (os.path.exists(file_dir + 'genericmask.txt')==True):
        print "Loading existing generic mask"
        generic_mask_file = file_dir + 'genericmask.txt'
        coords_generic = np.loadtxt(generic_mask_file)

    #Load Artifact Mask
    coords_artifact=[]
    if (os.path.exists(file_dir + 'artifactmask.txt')==True):
        print "Loading existing artifact mask"
        artifact_mask_file = file_dir + 'artifactmask.txt'
        coords_artifact = np.loadtxt(artifact_mask_file)

    #Load Bregma Mask
    coords_bregma=[]
    if (os.path.exists(file_dir + 'bregmamask.txt')==True):
        print "Loading existing bregma mask"
        bregma_mask_file = file_dir + 'bregmamask.txt'
        coords_bregma = np.loadtxt(bregma_mask_file)

    #PROCESS IMAGE TIMES FIRST
    #Test only on images in first X secs
    time_length = 50
    temp0 = np.where(np.logical_and(img_times>=img_times[0], img_times<=img_times[0]+time_length))[0]
    img_times = img_times[temp0]

    print "No. images: ", len(img_times)
    print "First and last img time: ", img_times[0], img_times[-1]
    #print len(img_times)

    selected_images = np.array(images_raw[temp0[:-1]],dtype=np.float32) #Remoe last image as it is not being use

    #remove baseline for [Ca]; but not VSD?
    if True:
        print "Removing baseline over all data"
        baseline = np.mean(selected_images, axis=0)
        selected_images = (selected_images - baseline)/baseline
    
    #PROCESS ephys data;
    sim_dir = file_dir
    sorted_file = file_name
    
    low_pass = False
    raw = True
    mua_load = False
    
    if low_pass:
        f = open(sim_dir + sorted_file+'/' + sorted_file + '_lp.tsf', "rb")
        tsf = Tsf_file(f, sim_dir)  #Auto load tsf file attributes: n_electrodes, ec_traces, SampleFreqeuncy and others
        tsf.sim_dir = sim_dir
        tsf.tsf_name = sim_dir + sorted_file+'/' + sorted_file + '.tsf'
        tsf.fname = sim_dir + sorted_file+'/' + sorted_file + '.tsf'
        f.close()
    
    if raw:
        f = open(sim_dir + sorted_file+'/' + sorted_file + '_raw.tsf', "rb")
        tsf = Tsf_file(f, sim_dir)  #Auto load tsf file attributes: n_electrodes, ec_traces, SampleFreqeuncy and others
        tsf.sim_dir = sim_dir
        tsf.tsf_name = sim_dir + sorted_file+'/' + sorted_file + '.tsf'
        tsf.fname = sim_dir + sorted_file+'/' + sorted_file + '.tsf'
        f.close()
    
    if mua_load:
        #Load Sorted data
        work_dir = sim_dir + file_name + "/"
        file_name = sorted_file + '_hp'
        ptcs_flag = 0
        Sort = Loadptcs(file_name, work_dir, ptcs_flag, save_timestamps=False)
        Sort.name=file_name
        Sort.filename=file_name
        Sort.directory=work_dir
       
        
    print "Sample freq: ", tsf.SampleFrequency
    print "Loaded .tsf for xcorrelation"
    
    #Load only ephys data within imaging window.
    #Find time indexes in ephys data that fall within imaging window - temp1

    if True:  #If using ephys data
        ephys_times = np.arange(0,float(len(tsf.ec_traces[0])),1.E3/float(tsf.SampleFrequency))*1.E-3
        temp1 = np.where(np.logical_and(ephys_times>=img_times[0], ephys_times<=img_times[-1]))[0]
        print "Ephys time indexes: ", np.float32(temp1)/tsf.SampleFrequency

        print "Computing split_array used to chunk the trigger data"
        split_array = []
        for i in range(1,len(img_times)-1,1):
            split_array.append(int(len_frame*i*tsf.SampleFrequency))

    bands = [[0.1, 4.0], [4.0, 8.0], [8.0, 12.0], [12.0, 25.0], [25.0, 100.0], [500.0, 5000.0]]
    bands_names = ['Delta (0.1-4.0)', 'Theta (4.0-8.0)', 'Alpha (8.0-12.0)', 'Beta (12.0-25.0)', 'Gamma (25.0-150.0)', 'High (500-5000)', 'MUA (Sorted)']
    
    counter = 0
    top_row = True
    electrodes = np.arange(0,16,1)
    electrodes = [0,5,10,15]
    n_electrodes = len(electrodes)
    
    for q in electrodes:
        for b in range(len(bands)):
            ephys_data = tsf.ec_traces[q][temp1].copy() * -1.0
            fs = 20000 
            lowcut, highcut = bands[b]
            ephys_data = butter_bandpass_filter(ephys_data, lowcut, highcut, fs, order = 2)

            print ""
            print ""
            print "Channel: ", q

            #print "Splitting ephys_data into chunks triggered on imaging times"
            ephys_split = np.split(ephys_data, split_array)
            ephys_split_mean = []
            for i in range(len(ephys_split)):
                ephys_split_mean.append(np.mean(ephys_split[i]))
                
            if False:  #Show LFP data overlayed with averaged chunks triggered on image times
                width = []
                for i in range(len(img_times)):
                    width.append(10)
                plt.bar(img_times, width, .01, color='black', alpha=0.65)
                plt.plot(ephys_times[temp1], ephys_data, color='blue')
                plt.bar(img_times[:-1], ephys_split_mean, .03, color='pink', alpha=0.45)
                plt.show()
            
            #print "Computing mean of ephys_split data"
            ephys_split_mean = np.array(ephys_split_mean)
            ephys_split_mean_max = np.max(np.abs(ephys_split_mean))
            ephys_split_mean = ephys_split_mean/ephys_split_mean_max
            ephys_split_mean = np.clip(ephys_split_mean,0,100) #remove negative values;
            
            #print "Computing image_lfp"  #MUST PARALLELIZE this OR FIGURE OUT HOW TO DO IT IN NUMPY ARRAYS
            
            image_lfp = np.einsum('m,mdr->mdr',ephys_split_mean, selected_images)
            image_lfp = np.mean(image_lfp, axis=0)

            #print "Plotting image_lfp"

            ax = plt.subplot(n_electrodes, 7, counter+1)
            ax.get_xaxis().set_visible(False)
            ax.set_yticklabels([])
            if top_row:
                plt.title(bands_names[counter])
                
            #Ablate generic map
            min_pixel = np.min(image_lfp)
            for j in range(len(coords_generic)):
                image_lfp[min(n_pixels,int(coords_generic[j][0]))][min(n_pixels,int(coords_generic[j][1]))]=min_pixel
            for j in range(len(coords_artifact)):
                image_lfp[min(n_pixels,int(coords_artifact[j][0]))][min(n_pixels,int(coords_artifact[j][1]))]=min_pixel
            for j in range(len(coords_bregma)):
                for k in range(n_pixels):
                    for l in range(7):
                        image_lfp[k][min(n_pixels-1,int(coords_bregma[1])-3+l)]=min_pixel
                        
            #for j in range(len(coords_bregma)):
            #    image_lfp[min(n_pixels,int(coords_bregma[j][0]))][min(n_pixels,int(coords_bregma[j][1]))]=min_pixel
            #    
            plt.imshow(image_lfp, origin='lower')
            if counter%7==0:
                plt.ylabel("Ch: "+str(q))
            plt.xlim(0,n_pixels-1)
            plt.ylim(n_pixels-1, 0)
            counter+=1

        #PLOT MUA AT END OF TURD
        print "Plotting MUA"
        ax = plt.subplot(n_electrodes, 7, counter+1)
        plt.xlim(0,n_pixels-1)
        plt.ylim(n_pixels-1, 0)
        ax.get_xaxis().set_visible(False)
        ax.set_yticklabels([])
        
        Sort_temp = Object_temp()
        Sort_temp.units = []
        Sort_temp.samplerate = Sort.samplerate
        for i in range(len(Sort.units)):
            if Sort.maxchan[i]==q:
                Sort_temp.units.append(Sort.units[i])
        
        #print len(Sort_temp.units)
        if len(Sort_temp.units)>0: 
            mua = MUA_compute(Sort_temp, len_frame)  #2D Array containing times and values of mua 
        
            #Compute indexes for mua
            mua_times = mua[0]
            temp2 = np.where(np.logical_and(mua_times>=img_times[0], mua_times<=img_times[-1]))[0]
            #print "Mua time indexes: ", ephys_times[temp1]#/len_frame

            mua_data = mua[1][temp2] #Load mua data 
            mua_split_array = []
            for i in range(1,len(img_times)-1,1):
                mua_split_array.append(int(i))
            
            mua_split = np.split(mua_data, mua_split_array)
            mua_split_mean = []
            for i in range(len(mua_split)):
                mua_split_mean.append(np.mean(mua_split[i]))

            if False:  #Show LFP data overlayed with averaged chunks triggered on image times
                plt.bar(img_times, width, .01, color='black', alpha=0.65)
                plt.plot(mua_times[temp2], mua_data, color='blue')
                plt.bar(img_times[:-1], mua_split_mean, .03, color='pink', alpha=0.45)
                plt.show()
            
            mua_split_mean = np.array(mua_split_mean)
            mua_split_mean_max = np.max(np.abs(mua_split_mean))
            if mua_split_mean_max>0: mua_split_mean = mua_split_mean/mua_split_mean_max
            mua_split_mean = np.clip(mua_split_mean,0,100) #remove negative values;
            
            image_lfp = np.einsum('m,mdr->mdr',mua_split_mean, selected_images)
            image_lfp = np.mean(image_lfp, axis=0)
            
            #Ablate generic and artifact maps
            min_pixel = np.min(image_lfp)
            for j in range(len(coords_generic)):
                image_lfp[min(n_pixels,int(coords_generic[j][0]))][min(n_pixels,int(coords_generic[j][1]))]=min_pixel
            for j in range(len(coords_artifact)):
                image_lfp[min(n_pixels,int(coords_artifact[j][0]))][min(n_pixels,int(coords_artifact[j][1]))]=min_pixel
            for j in range(len(coords_bregma)):
                for k in range(n_pixels):
                    for l in range(7):
                        image_lfp[k][min(n_pixels-1,int(coords_bregma[1])-3+l)]=min_pixel
            
            plt.imshow(image_lfp, origin='lower')
        else:
            plt.imshow(image_lfp*0.0, origin='lower')

        if top_row:
            plt.title(bands_names[counter])

        counter+=1
        top_row=False
    
    #Plot EEG DATA
    
    
    plt.suptitle(file_name, fontsize = 20)
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()

    #plt.tight_layout()
    plt.show()
    quit()
    


def Multitaper_specgram_allfreqs(time_series):

    s = time_series
    
    #print "Length of recording: ", len(s)
    
    #******************************************************
           
    #Multi taper function parameters
    nfft = 5000
    shift = 10
    Np = 1000
    k = 6
    tm = 6.0
    
    #Computing multi taper trf_specgram
    spec = tfr_spec(s, nfft, shift, Np, k, tm)

    #Martin changes
    zis = np.where(spec == 0.0) # row and column indices where P has zero power
    if len(zis[0]) > 0: # at least one hit
        spec[zis] = np.finfo(np.float64).max # temporarily replace zeros with max float
        minnzval = spec.min() # get minimum nonzero value
        spec[zis] = minnzval # replace with min nonzero values
    spec = 10. * np.log10(spec) # convert power to dB wrt 1 mV^2?

    p0=-40
    p1=None
    if p0 != None:
        spec[spec < p0] = p0
    if p1 != None:
        spec[spec > p1] = p1

    return spec, nfft
    
def Compute_lpf_static_maps(file_dir, file_name, images_raw, img_start, img_end, len_frame, img_rate, n_pixels, n_frames, img_times):

    #Load Generic Mask
    coords_generic=[]
    if (os.path.exists(file_dir + 'genericmask.txt')==True):
        print "Loading existing generic mask"
        generic_mask_file = file_dir + 'genericmask.txt'
        coords_generic = np.loadtxt(generic_mask_file)

    #Load Artifact Mask
    coords_artifact=[]
    if (os.path.exists(file_dir + 'artifactmask.txt')==True):
        print "Loading existing artifact mask"
        artifact_mask_file = file_dir + 'artifactmask.txt'
        coords_artifact = np.loadtxt(artifact_mask_file)

    #Load Bregma Mask
    coords_bregma=[]
    if (os.path.exists(file_dir + 'bregmamask.txt')==True):
        print "Loading existing bregma mask"
        bregma_mask_file = file_dir + 'bregmamask.txt'
        coords_bregma = np.loadtxt(bregma_mask_file)

    #PROCESS IMAGE TIMES FIRST
    #Test only on images in first X secs
    time_length = 120    #No. of secs of analysis
    temp0 = np.where(np.logical_and(img_times>=img_times[0], img_times<=img_times[0]+time_length))[0]
    img_times = img_times[temp0]

    print "No. images: ", len(img_times)
    print "First and last img time: ", img_times[0], img_times[-1]
    #print len(img_times)

    #selected_images = np.array(images_raw[temp0[:-1]],dtype=np.float32) #Remove last image as it is not being use
    selected_images = np.array(images_raw[temp0],dtype=np.float32) #Remove last image as it is not being use

    #remove baseline for [Ca]; but not VSD?
    if True:
        print "Removing baseline over all data"
        baseline = np.mean(selected_images, axis=0)
        selected_images = (selected_images - baseline)/baseline
    
    #PROCESS ephys data;
    sim_dir = file_dir
    sorted_file = file_name
    
    low_pass = False
    raw = True
    mua_load = False
    
    if low_pass:
        f = open(sim_dir + sorted_file+'/' + sorted_file + '_lp.tsf', "rb")
        tsf = Tsf_file(f, sim_dir)  #Auto load tsf file attributes: n_electrodes, ec_traces, SampleFreqeuncy and others
        tsf.sim_dir = sim_dir
        tsf.tsf_name = sim_dir + sorted_file+'/' + sorted_file + '.tsf'
        tsf.fname = sim_dir + sorted_file+'/' + sorted_file + '.tsf'
        f.close()
    
    if raw:
        f = open(sim_dir + sorted_file+'/' + sorted_file + '_raw.tsf', "rb")
        tsf = Tsf_file(f, sim_dir)  #Auto load tsf file attributes: n_electrodes, ec_traces, SampleFreqeuncy and others
        tsf.sim_dir = sim_dir
        tsf.tsf_name = sim_dir + sorted_file+'/' + sorted_file + '.tsf'
        tsf.fname = sim_dir + sorted_file+'/' + sorted_file + '.tsf'
        f.close()
    
    if mua_load:
        #Load Sorted data
        work_dir = sim_dir + file_name + "/"
        file_name = sorted_file + '_hp'
        ptcs_flag = 0
        Sort = Loadptcs(file_name, work_dir, ptcs_flag, save_timestamps=False)
        Sort.name=file_name
        Sort.filename=file_name
        Sort.directory=work_dir
       
        
    print "Sample freq: ", tsf.SampleFrequency
    print "Loaded .tsf for xcorrelation"
    
    #Load only ephys data within imaging window.
    #Find time indexes in ephys data that fall within imaging window - temp1

    if True:  #If using ephys data
        print "downsampling ephys to 1khz...", len(tsf.ec_traces[0])
        ephys_times = np.arange(0,float(img_times[-1])*1E3,1.E3/float(tsf.SampleFrequency))*1.E-3
        print "finding matching imaging frames for ephys times..."
        temp1 = np.where(np.logical_and(ephys_times>=img_times[0], ephys_times<=img_times[-1]))[0]
        print "Ephys time indexes: ", np.float32(temp1)/tsf.SampleFrequency

        #computing split_array used to chunk the ephys data: used to obtain average LFP during single img frame; timesteps of original data
        split_array = []
        for i in range(0,len(img_times)-1,1):
            split_array.append(int(len_frame*i*tsf.SampleFrequency))

    bands = [[0.1, 4.0], [4.0, 8.0], [8.0, 12.0], [12.0, 25.0], [25.0, 100.0], [500.0, 5000.0]]
    #bands = [[0.1,5000.]]
    bands_names = ['Delta (0.1-4.0)', 'Theta (4.0-8.0)', 'Alpha (8.0-12.0)', 'Beta (12.0-25.0)', 'Gamma (25.0-150.0)', 'High (500-5000)']
    #bands_names = ["all"]
    bands_out = ['delta', 'theta', 'alpha', 'beta', 'gamma', 'high']
    #bands_out = ['all']
    
    counter = 0
    top_row = True
    #electrodes = np.arange(0,16,1)
    electrodes = np.arange(0,16,1)
    #electrodes = [0,5,10,15]
    
    n_electrodes = len(electrodes)
    
    lfp_correlation_method = True
    lfp_power_method = False    #Requires split array +1 and selected_images - 1... not clear why.
    
    for q in electrodes:
        if lfp_correlation_method:
            for b in range(len(bands)):
                print "Channel: ", q, " band: ", bands_names[b]
                
                #Use LFP Correlation method: average lfp signal during each frame and mutiply by frame
                ephys_data = tsf.ec_traces[q][temp1].copy() * 1.0       #Can look at only postivie or neative data
                fs = 20000 
                lowcut, highcut = bands[b]
                ephys_data = butter_bandpass_filter(ephys_data, lowcut, highcut, fs, order = 2)

                ephys_data = np.clip(ephys_data, 0, 1E10) #Look only at positive power ;
                #ephys_data = np.clip(ephys_data, -1E10, 0) #Look only at negative power;
                
                #print splitting ephys_data into chunks triggered on imaging times: to get average value of LFP during single img frame
                ephys_split = np.split(ephys_data, split_array) #
                
                ephys_split_mean = []
                for i in range(len(ephys_split)):
                    ephys_split_mean.append(np.mean(ephys_split[i]))

                    
                ##Visualize img frame-to-ephys matching: LFP data overlayed with averaged chunks triggered on image times
                #if True:
                    #width = []
                    #for i in range(len(img_times)):
                        #width.append(10)
                    #plt.bar(img_times, width, .01, color='black', alpha=0.65)
                    #plt.plot(ephys_times[temp1], ephys_data, color='blue')
                    #print len(img_times), ephys_split_mean.shape
                    #plt.bar(img_times, ephys_split_mean, .03, color='pink', alpha=0.45)
                    #plt.show()
                
                #print "Computing mean of ephys_split data"
                ephys_split_mean = np.clip(ephys_split_mean,0,1000000) #remove negative values;
                ephys_split_mean = np.float32(ephys_split_mean)
                ephys_split_mean_max = np.nanmax(ephys_split_mean)     #Normalize to largest LFP value (negative or positive)
                ephys_split_mean = ephys_split_mean/ephys_split_mean_max    
                #ephys_split_mean = np.clip(ephys_split_mean,0,1) #remove negative values;
                ephys_split_mean = np.nan_to_num(ephys_split_mean)

                #Compute lfp average triggered images
                image_lfp = np.einsum('m,mdr->mdr',ephys_split_mean, selected_images)
                image_lfp = np.mean(image_lfp, axis=0)

                np.save(file_dir + file_name+'/'+file_name+'_lfpmap_band_'+bands_out[b]+'_channel_'+str(q).zfill(2), image_lfp)

        if lfp_power_method:
            ephys_data = tsf.ec_traces[q][temp1].copy()     #This selects only part of recording
            fs = 20000 
            
            #lowcut = 0.1
            #highcut = 110.
            #ephys_temp = butter_bandpass_filter(ephys_data, lowcut, highcut, fs, order = 2)
            print "Computing time-frequency reassignment specgram channel: ", q
            
            tfr_file = file_dir + file_name+'/'+file_name+'_tfr_channel_'+str(q).zfill(2)
            if os.path.exists(tfr_file+'.npz'): 
                data = np.load(tfr_file+'.npz')
                mt_specgram = data['mt_specgram']
                nfft = data['nfft']
            else: 
                mt_specgram, nfft = Multitaper_specgram_allfreqs(ephys_data)
                np.savez(tfr_file, mt_specgram=mt_specgram, nfft=nfft)
            
            #plt.imshow(mt_specgram[::-1][:10000,:10000], origin='upper', aspect='auto', interpolation='none')
            #plt.show()
            
            #SampleFrequency = tsf.samplerate
            #mt_specgram, extent = Compute_specgram_signal(ephys_data, SampleFrequency):
            

            for b in range(len(bands)):
                f0, f1 = bands[b]
                lo = int(float(nfft)/1000. * f0)
                hi = int(float(nfft)/1000. * f1)
                mt_specgram_temp = mt_specgram[lo:hi][::-1] #Take frequency band slice only; also invert the data

                #print f0, f1
                #plt.imshow(mt_specgram_temp, extent = [0,f1, len(mt_specgram_temp),0], origin='upper', aspect='auto', interpolation='sinc')
                #plt.show()
            
                spec_split = np.array(split_array)*len(mt_specgram_temp[0])/split_array[-1] #Normalize splitting array to length of spectrogram
                #print spec_split
                spec_split_mean = []
                for ss in range(len(spec_split)-1):
                    spec_split_mean.append(np.mean(mt_specgram_temp[:,spec_split[ss]:spec_split[ss+1]]))

                #print "Computing mean of ephys_split data"
                spec_split_mean = np.array(spec_split_mean)
                spec_split_mean_max = np.max(np.abs(spec_split_mean))     #Normalize to largest LFP value (negative or positive); THIS CAN"T BE NEGATIVE!
                spec_split_mean_min = np.min(np.abs(spec_split_mean))     #Normalize to largest LFP value (negative or positive); THIS CAN"T BE NEGATIVE!
                spec_split_mean = (spec_split_mean-spec_split_mean_min)/(spec_split_mean_max-spec_split_mean_min)
                #plt.plot(spec_split_mean)
                #plt.show()
                #quit()
                #spec_split_mean = np.clip(spec_split_mean,0,1) #remove negative values; NOT REQUIRED

                #Compute lfp average triggered images
                print spec_split_mean.shape, selected_images.shape
                image_lfp = np.einsum('m,mdr->mdr', spec_split_mean, selected_images)
                image_lfp = np.mean(image_lfp, axis=0)

                np.save(file_dir + file_name+'/'+file_name+'_powermap_band_'+bands_out[b]+'_channel_'+str(q).zfill(2), image_lfp)

            #print "Plotting image_lfp"
            ax = plt.subplot(n_electrodes, 7, counter+1)
            ax.get_xaxis().set_visible(False)
            ax.set_yticklabels([])
            if top_row:
                plt.title(bands_names[counter])

            #Ablate generic map
            min_pixel = np.min(image_lfp)
            for j in range(len(coords_generic)):
                image_lfp[min(n_pixels,int(coords_generic[j][0]))][min(n_pixels,int(coords_generic[j][1]))]=min_pixel
            for j in range(len(coords_artifact)):
                image_lfp[min(n_pixels,int(coords_artifact[j][0]))][min(n_pixels,int(coords_artifact[j][1]))]=min_pixel
            for j in range(len(coords_bregma)):
                for k in range(n_pixels):
                    for l in range(7):
                        image_lfp[k][min(n_pixels-1,int(coords_bregma[1])-3+l)]=min_pixel
                        
            #for j in range(len(coords_bregma)):
            #    image_lfp[min(n_pixels,int(coords_bregma[j][0]))][min(n_pixels,int(coords_bregma[j][1]))]=min_pixel
            #    
            plt.imshow(image_lfp, origin='lower')
            if counter%7==0:
                plt.ylabel("Ch: "+str(q))
            plt.xlim(0,n_pixels-1)
            plt.ylim(n_pixels-1, 0)
            counter+=1


        counter+=1
        top_row=False
   
    #plt.suptitle(file_name, fontsize = 20)
    #figManager = plt.get_current_fig_manager()
    #figManager.window.showMaximized()

    #plt.show()
    #quit()
    

def Compute_specgram_signal(data, SampleFrequency, f0=0.1, f1=110, p0=-40):

    t0=None
    t1=None
    #f0=0.1
    #f1=110
    #p0         #clipping bottom of specgram
    p1=None
    chanis=-1
    width=2
    tres=.5
    cm=None
    colorbar=False
    title=True
    figsize=(20, 6.5)
    
    P0, P1 = None, None
    chanis = -1

    sampfreq=SampleFrequency #1KHZ LFP SAMPLE RATE for Nick's data; Otherwise full sample rates;

    #NFFT = intround(width * sampfreq)
    #NOVERLAP = intround(NFFT - tres * SAMPFREQ)

    length = len(data)

    ts = np.arange(0,len(data),1.0)/sampfreq

    if t0 == None:
        t0, t1 = ts[0], ts[-1] # full duration
    #if t1 == None:
    #    t1 = t0 + 10 # 10 sec window
    if width == None:
        width = uns['LFPWIDTH'] # sec
    if tres == None:
        tres = uns['LFPTRES'] # sec
    assert tres <= width

    NFFT = intround(width * sampfreq)
    noverlap = intround(NFFT - tres * sampfreq)

    t0i, t1i = ts.searchsorted((t0, t1))

    #data = filter.notch(data)[0] # remove 60 Hz mains noise

    print "Computing regular fft specgram"
    P, freqs, t = mpl.mlab.specgram(data/1e3, NFFT=NFFT, Fs=sampfreq, noverlap=noverlap)
    
    # convert t to time from start of acquisition:
    t += t0
    # keep only freqs between f0 and f1:
    if f0 == None:
        f0 = freqs[0]
    if f1 == None:
        f1 = freqs[-1]
    lo, hi = freqs.searchsorted([f0, f1])
    P, freqs = P[lo:hi], freqs[lo:hi]
    #print P
    
    # check for and replace zero power values (ostensibly due to gaps in recording)
    # before attempting to convert to dB:
    zis = np.where(P == 0.0) # row and column indices where P has zero power
    if len(zis[0]) > 0: # at least one hit
        P[zis] = np.finfo(np.float64).max # temporarily replace zeros with max float  #CAT: This can probably be unhacked using nanmax or masked arrays
        minnzval = P.min() # get minimum nonzero value
        P[zis] = minnzval # replace with min nonzero values
    P = 10. * np.log10(P) # convert power to dB wrt 1 mV^2?

    # for better visualization, clip power values to within (p0, p1) dB
    if p0 != None:
        P[P < p0] = p0
    if p1 != None:
        P[P > p1] = p1

    extent = ts[0], ts[-1], freqs[0], freqs[-1]

    return P[::-1], extent
    
    

def load_traces(stmtd, n_areas, window, area_names, maxmin, hemisphere, exp_dirs, sub_dirs, main_dir):
        
    cell_counter=0
    exp_counter =0
    n_interp_points = 600
    normalize=False

    for exp_dir in exp_dirs:
        print exp_dir
        
        sub_dir_counter=0
        temp_counter=0
        for sub_dir in sub_dirs[exp_counter]:
            print sub_dir

            #Load depth file (cortex vs. subcortical)
            depth_file = exp_dir+'/'+sub_dir+'/'+'depth.txt'
            with open(depth_file, "r") as f:
                depth_temp = f.read().splitlines()[0]

            #Load cortical state (awake vs. anesthetized)
            state_file = exp_dir+'/'+sub_dir+'/'+'state.txt'
            with open(state_file, "r") as f:
                state_temp = f.read().splitlines()[0]

            #Load rasters and channel of all cells in experiment
            work_dir = exp_dir + "/"+sub_dir+'/'
            file_name =  sub_dir+"_hp" #.ptcs"
            Sort = Loadptcs(file_name, work_dir, 1, save_timestamps=False) #Auto load flag for Nick's data
            Sort.name=file_name
            Sort.filename=file_name
            Sort.directory=work_dir
            
            #Computing or loading cell template from file
            width_template = 40
            if (os.path.exists(exp_dir + "/"+sub_dir+'/templates.npy')==False):
                tsf_name = exp_dir + "/"+sub_dir+'/'+sub_dir+'_hp.tsf'
                f = open(tsf_name, "rb")
                tsf = Tsf_file(f, exp_dir + "/"+sub_dir+'/')  #Auto load tsf file attributes: n_electrodes, ec_traces, SampleFreqeuncy and others
                tsf.sim_dir = exp_dir + "/"+sub_dir+'/'
                tsf.tsf_name = tsf_name
                f.close()

                raw_templates= []
                temp_ray = np.zeros(width_template*2,dtype=np.float32)
                for i in range(len(Sort.units)):
                    temp_ray=temp_ray*0.0
                    #print "no. spikes in cell: ", len(Sort.units[i])
                    for spike in Sort.units[i]:
                        if (spike>width_template) and (spike<(Sort.units[i][-1]-width_template)):  #Need to ensure spike not too close to beginning or end of recording
                            temp_ray+= tsf.vscale_HP*(tsf.ec_traces[Sort.maxchan[i]][spike-width_template:spike+width_template])
                    temp_ray = temp_ray/len(Sort.units[i])
                    raw_templates.append(temp_ray)
                np.save(exp_dir + "/"+sub_dir+'/templates', raw_templates)
                
            else:
                raw_templates = np.load(exp_dir + "/"+sub_dir+'/templates.npy')
            
            #List all files in exp directory
            file_dirs = os.listdir(exp_dir+'/'+sub_dir+'/')
            for file_ in file_dirs:
                #if ("time_course_data" in file_):# and ("all" in  file_):
                if ("time_course_data" in file_) and (area_names[0] in file_):
                    temp_file = file_.replace("time_course_data_"+sub_dir+"_all_3sec_window_unit",'')
                    unit = temp_file[0:2]
                    
                    #Use csv_file_name to load ptp and n_spikes variables
                    csv_file_name = glob.glob(exp_dir+'/'+sub_dir+'/'+sub_dir+'_all_maxmap_unit_'+unit+'*.npy')
                    #print exp_dir+'/'+sub_dir+'/'+sub_dir+'_all_maxmap_unit_'+unit+'*.npy'
                    #print csv_file_name
                    stmtd.ptp.append(int(csv_file_name[0][-20:-17]))
                    stmtd.n_spikes.append(int(csv_file_name[0][-9:-4]))

                    #Check to see if sufficient spikes in cell (first value in txt file)
                    #if ptp>ptp_cutoff:
                        
                    #Load time-courses for all cortical areas; data[1] = min value; data[2] = max value
                    data=[]
                    with open(exp_dir+'/'+sub_dir+'/'+file_, "r") as f:
                        temp_data = csv.reader(f)
                        for row in temp_data:
                            data.append(row)
                    
                    temp_counter+=1

                    #Interpolate time course as there are different imaging rates across the experiments
                    x = np.arange(0, len(data[2]), 1)
                    #print data[2]
                    #print data[1]
                    if (np.max(np.float32(data[2]))>-np.min(np.float32(data[1]))):
                        f = interpolate.interp1d(x, data[2])
                    else:
                        f = interpolate.interp1d(x, data[1])
                    xx = np.linspace(x[0],x[-1],n_interp_points)
                    data_max = f(xx)

                    stmtd.tc_data.append(np.float32(data_max)*1E2)  #Load max cortex pixel time course

                    #Ensure signal > than minDF_response limit:; 2 cases: (i) DF/F > 1%; or (ii) n_spikes>cutoff and 
                    max_DFF_response = max(max(data_max),abs(min(data_max)))
                    stmtd.DF.append(max_DFF_response)

                    #Normalize channel location 
                    stmtd.channel_left.append(float(Sort.maxchan[int(unit)]-min(Sort.maxchan))/(max(Sort.maxchan)-min(Sort.maxchan)))
                    stmtd.raster_left.append(Sort.units[int(unit)])
                    stmtd.template_left.append(raw_templates[int(unit)])

                    #Save index of trace area
                    stmtd.area_index_left.append(area_names[0])
                    
                    #Save name of trace if "max" or "min"
                    stmtd.area_maxormin.append('max')
                    
                    #Save hemisphere of trace
                    stmtd.area_hemisphere.append('left')

                    #Save depth of cell
                    stmtd.depth_left.append(depth_temp)
                    
                    #Save cortical state of cell
                    stmtd.state_left.append(state_temp)
                    
                    #Mouse counter: label the trace with unique mouse label
                    stmtd.mouse_index.append(exp_counter)
                    
                    stmtd.file_index.append(file_)
                    stmtd.dir_index.append(exp_dir+sub_dir)

                    cell_counter+=1
            sub_dir_counter+=1
        exp_counter+=1
        print "# cells in experiment: ", temp_counter
    print "Total # cells: ", cell_counter
    print ""
    print ""
    #print sys.getsizeof(stmtd)

def recombine_data():
    ''' Function to recombine parcelated matrices from  sta_clustering code
    '''
    #recombine data
    n_rows = int(len(data[0])/size)
    n_cols = int(len(data[0])/size)
    data_blank = np.zeros((len(data), size,size))-0.005
    
    data_array = []
    for k in range(0,n_rows,1):
        temp_stack = data_split[k*n_rows]
        for p in range(1, n_cols, 1):
            if ((k+p) % 2)==0:
                #temp_stack = np.dstack((temp_stack, data_blank))       #Parcelation
                temp_stack = np.ma.dstack((temp_stack, data_split[k*n_rows+p]))
            else:
                temp_stack = np.ma.dstack((temp_stack, data_split[k*n_rows+p]))

        if len(data_array)==0:
            data_array = np.ma.array(temp_stack)
        else:
            data_array = np.ma.hstack((data_array, temp_stack))
    
    data = np.ma.array(data_array)
    print data.shape

def animate_data(data):
    '''Make movies from data matrix
    '''
    v_max=np.nanmax(data)
    v_min=np.nanmin(data)
    print v_max, v_min
    
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=5, metadata=dict(artist='Me'), bitrate=1800)

    fig = plt.figure() # make figure

    # make axesimage object
    im = plt.imshow(data[0], cmap=plt.get_cmap('jet'), vmin=v_min, vmax=v_max, interpolation='none')#, vmin=0, vmax=v_max)
    # function to update figure
    def updatefig(j):
        # set the data in the axesimage object
        im.set_array(data[j])
        plt.title("Frame: "+str(j)+"\n"+str(round(float(j)/150,2))+"sec")
        # return the artists set
        return im,
    # kick off the animation
    ani = animation.FuncAnimation(fig, updatefig, frames=range(len(data)), interval=100, blit=False, repeat=True)

    if False:
        ani.save(fname+'.mp4', writer=writer)
    
    plt.show()

def mask_data(data, main_dir, n_pixels):
            
    #Load General mask (removes background)
    generic_mask_file = []
    generic_mask_file = main_dir + 'genericmask.txt'
    if (os.path.exists(generic_mask_file)==False):
        generic_coords = Define_generic_mask(data, main_dir)
    else:
        generic_coords = np.loadtxt(generic_mask_file)
        
    generic_mask_indexes=np.zeros((n_pixels,n_pixels))
    for i in range(len(generic_coords)):
        generic_mask_indexes[generic_coords[i][0]][generic_coords[i][1]] = True
    
    temp_array = np.ma.array(np.zeros((len(data),n_pixels,n_pixels),dtype=np.float32), mask=True)
    #Mask all frames; NB: PROBABLY FASTER METHOD
    for i in range(0, len(data),1):
        temp_array[i] = np.ma.masked_array(data[i], mask=generic_mask_indexes)
    
    return temp_array
    
def parcelation_clustering():
    area_names = ['hindlimb', 'forelimb', 'whisker','retrosplenial','visual'] 
    #area_names = area_names[::-1]
    
    dir_counter = 0
    for file_dir in file_dirs:
        
        for file_name in file_names[dir_counter]:

            #Load location of recording:
            depth = np.loadtxt(file_dir+file_name+'/depth.txt', dtype=str)
            print file_name, " ", depth

            #Load anesthetic state:
            state = np.loadtxt(file_dir+file_name+'/state.txt', dtype=str)
            #if state != state_match: continue

            #Load units from .csv file name; NB: This may lead to duplication as resaved .csv's for Dongsheng
            files = os.listdir(file_dir+file_name)
            temp_names = []
            for file_ in files:
                if ("unit_" in file_) and ('.csv' in file_):
                    temp_names.append(file_)

            #Save individual unit names, channels, ptps
            units = []
            channels = []
            ptps = []
            for file_ in temp_names:
                units.append(int(file_[5:7]))
                channels.append(int(file_[16:18]))
                ptps.append(int(file_[23:26]))

            for i in range(len(units)):
                #print "Unit: ", units[i]
                #Load .csv file to get # spikes
                csv_file = glob.glob(file_dir + file_name + "/unit_"+str(units[i]).zfill(2) + "*.csv")[0]
                spikes = np.loadtxt(csv_file, dtype=np.float32)
                
                if (ptps[i] < 35.0) or (len(spikes)<128): continue

                #************** Load time course
                #print file_dir+file_name+'/img_avg_'+file_name+"_unit"+str(units[i])+"*.npy"
                
                parcel_file = file_dir+file_name+"/"+file_name+'_parcelation_'+str(25)+"_unit_"+str(units[i]).zfill(2)
                if(os.path.exists(parcel_file)==False):
                    
                    unit_file = glob.glob(file_dir+file_name+'/img_avg_'+file_name+"_unit"+str(units[i]).zfill(2)+"_ch"+str(channels[i])
                    +"_all_3sec*spikes.npy")
                    if len(unit_file)==0: continue
                    else: unit_file = unit_file[0]
                    
                    
                    #print "Loading ...", unit_file
                    data = np.load(unit_file)
                    print data.shape

                    #**************** Mask Data *****************
                    print "Masking data..."
                    data = mask_data(data, main_dir, 256)
                    print data.shape
                    
                    #plt.imshow(data[0])
                    #plt.show()

                    #****************Split data into parcels
                    print "Splitting data..."
                    size = 25
                    n_rows = int(len(data[0])/size)
                    n_cols = int(len(data[0])/size)
                    
                    data_split = np.ma.array(np.zeros((n_rows*n_cols,len(data),size,size), dtype=np.float32), mask=True)
                    for p in range(0, n_rows, 1):
                        for j in range(0, n_cols, 1):
                            data_split[p*n_rows+j] = data[:,p*size:(p+1)*size,j*size:(j+1)*size]
                    
                    print data_split.shape

                    #******RECOMBINE DATA
                    if False: recombine_data()

                    #************ Animation display
                    if False:  animate_data(data)
                    
                    #************ Average activity 
                    print "Averaging data..."
                    data_ave = np.ma.average(np.ma.average(data_split, axis=2),axis=2)
                    
                    mask_indexes = []
                    for k in range(len(data_ave)):
                        if data_ave[k].mask.all(): mask_indexes.append(k)
                    
                    #print mask_indexes
                    #print len(mask_indexes)
                    #quit()
                                        
                    #data_ave = np.array(data_ave)
                    #data_ave = data_split
                    #print data_ave.shape
                    print ""
                    print ""
                    
                    data_ave.dump(parcel_file)

                    #quit()
                    
        dir_counter+=1
        
        
    quit()

def Meanshift(data):
    from sklearn.cluster import MeanShift, estimate_bandwidth
    from sklearn.datasets.samples_generator import make_blobs
    
    # The following bandwidth can be automatically detected using
    bandwidth = estimate_bandwidth(data, quantile=0.2, n_samples=1000)

    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    ms.fit(data)
    labels = ms.labels_
    cluster_centers = ms.cluster_centers_

    labels_unique = np.unique(labels)
    n_clusters_ = len(labels_unique)

    print("number of estimated clusters : %d" % n_clusters_)

    return labels

def plot_dendrograms():
    # Compute and plot first dendrogram.
    D_cortex = np.corrcoef(total_cortex_all)
    Y = sch.linkage(D_cortex, method='centroid')
    Z1_cortex = sch.dendrogram(Y, orientation='bottom')
    Y = sch.linkage(D_cortex, method='single')
    Z2_cortex = sch.dendrogram(Y)

    #Cortex plot
    ax1 = plt.subplot(121)
    img_out = np.corrcoef(total_cortex_all)
    plt.imshow(img_out, origin='lower', interpolation='none', cmap=plt.get_cmap('jet'))
    #n_cells = 0
    #for i in range(len(mouse_cortex_counter)):
        #n_cells += mouse_cortex_counter[i]
        #plt.plot([-0.5,len(total_cortex_all)-0.5], [n_cells-0.5,n_cells-0.5],'r--', color='black', linewidth=3)
    plt.ylim(-0.5,-0.5+len(D_cortex))
    plt.xlim(-0.5,-0.5+len(D_cortex))
    plt.title(depth_match + " " + state_match + "\n (in order of experiments)")

    #Distance clusters
    axmatrix = plt.subplot(122)
    idx1 = Z1_cortex['leaves']
    idx2 = Z2_cortex['leaves']
    D = D_cortex[idx1,:]
    D = D[:,idx2]
    im = axmatrix.matshow(D, origin='lower',interpolation='none', cmap=plt.get_cmap('jet'))
    plt.title(depth_match + " " + state_match + "\n(in cluster order)")
    #axmatrix.set_xticks(idx2)
    #axmatrix.set_yticks(idx1)

    xx = np.linspace(0,len(D_cortex)-1,len(D_cortex))
    plt.xticks(xx, idx2)
    plt.xticks(rotation=90)
    plt.yticks(xx, idx1)

    plt.ylim(-.5,-.5+len(D_cortex))
    plt.xlim(-.5,-.5+len(D_cortex))

    plt.show()
    
    quit()

def parcelation_1D(main_dir):
        
    #state_match = "anesthetized"
    state_match = "awake"
    
    #depth_match = "cortex"
    depth_match = "subcortical"
    
    dir_counter = 0
    vectors = []
    
    out_file = main_dir+"vectors_"+state_match+"_"+depth_match
    if(os.path.exists(out_file+'.npy')==False):
        for file_dir in file_dirs:
            for file_name in file_names[dir_counter]:

                #Load location of recording:
                depth = np.loadtxt(file_dir+file_name+'/depth.txt', dtype=str)
                print file_name, " ", depth
                if depth != depth_match: continue
                
                #Load anesthetic state:
                state = np.loadtxt(file_dir+file_name+'/state.txt', dtype=str)
                if state != state_match: continue

                #Load units from .csv file name; NB: This may lead to duplication as resaved .csv's for Dongsheng
                files = os.listdir(file_dir+file_name)
                temp_names = []
                for file_ in files:
                    if ("unit_" in file_) and ('.csv' in file_):
                        temp_names.append(file_)

                #Save individual unit names, channels, ptps
                units = []
                channels = []
                ptps = []
                for file_ in temp_names:
                    units.append(int(file_[5:7]))
                    channels.append(int(file_[16:18]))
                    ptps.append(int(file_[23:26]))

                for i in range(len(units)):
                    #Load .csv file to get # spikes
                    csv_file = glob.glob(file_dir + file_name + "/unit_"+str(units[i]).zfill(2) + "*.csv")[0]
                    spikes = np.loadtxt(csv_file, dtype=np.float32)
                    
                    if (ptps[i] < 35.0) or (len(spikes)<128): continue

                    
                    #************** Load Average activity data
                    unit_file = glob.glob(file_dir+file_name+'/'+file_name+"_parcelation_25_unit_"+str(units[i]).zfill(2))
                    if len(unit_file)==0: continue
                    else: unit_file = unit_file[0]
                    
                    data = np.load(unit_file)
                    
                    cell_vector=[]
                    for k in range(len(data)):
                        if (data[k].mask.all()==False): 
                            cell_vector.extend(data[k])
                    
                    x = np.linspace(0,len(cell_vector), len(cell_vector))
                    y = cell_vector
                    f = interpolate.interp1d(x,y)
                    x2 = np.linspace(0, len(cell_vector), 15000)
                    unrolled = f(x2)
                        
                    vectors.append(unrolled)
                    
                    #plt.ylim(-0.01, 0.01)
                    #plt.plot(unrolled)
                    #plt.show()
            
            dir_counter+=1 #Go to next mouse
        
        np.save(out_file, vectors)
    else:
        vectors = np.load(out_file+'.npy')


    #short_vector = []
    #marker = 15000./77.
    #for a in range(len(vectors)):
        #temp = []
        #for r in range(77):
            #temp.extend(vectors[a][(r+1/3.)*marker:(r+1)*marker])
        #short_vector.append(temp)
        
    #vector = short_vector

    #total_cortex_all = vectors

    return vectors


def plot_3D(labels):
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    n = 100
    c='r'
    m='^'
       

    clrs = []
    xs=[]
    ys=[]
    zs=[]

    cmap = mpl.cm.get_cmap('jet')

    n_clusters = 20 #n_clusters_

    clrs_rand = []
    for i in range(n_clusters):
        clrs_rand.append(np.random.rand(3,))
        
    for i in range(len(labels)):
        if labels[i]<n_clusters:
            
            clrs.append(colors[labels[i]])

            xs.append(data[i][0])
            ys.append(data[i][1])
            zs.append(data[i][2])
        
    ax.scatter(xs, ys, zs, c=clrs, marker=m, s = 100)
    plt.title(depth_match + "  " + state_match)
    plt.show()

def roi_clustering(data, file_dirs, file_names, main_dir, state_match, depth_match):
    
    if (os.path.exists(main_dir + 'total_'+depth_match+'_all_'+state_match+'.npy')==False):
    #if True:
        area_names = ['hindlimb', 'forelimb', 'barrel','retrosplenial','visual', 'motor', 'pta', 'acc'] 
        #area_names = area_names[::-1]
        
        left_tag_cortex_excitatory = []
        left_tag_cortex_inhibitory = []
        left_tag_subcortical_excitatory = []
        left_tag_subcortical_inhibitory = []
        
        
        #Make total roi arrays for each ROI and activation/inhibition types
        total_cortex_roi_left = []
        total_subcortical_roi_left = []
        for i in range(len(area_names)):
            total_cortex_roi_left.append([])
            total_subcortical_roi_left.append([])
            for j in range(2):
                total_cortex_roi_left[i].append([])
                total_subcortical_roi_left[i].append([])
        
        total_cortex_roi_left_cumulative= [[],[]]
        total_subcortical_roi_left_cumulative= [[],[]]
        
        total_cortex_all=[]
        total_cortex_all_files=[]
        total_cortex_all_ch=[]
        
        total_subcortical_all=[]
        total_subcortical_all_files=[]
        total_subcortical_all_ch=[]
        
        mouse_cortex_counter = []  #Keeps track of # of cells from each mouse
        mouse_subcortical_counter = []  #Keeps track of # of cells from each mouse
        
        dir_counter = 0
        for file_dir in file_dirs:
            
            cortex_unit_counter = 0 #Tracks # units in each mouse
            subcortical_unit_counter = 0 #Tracks # units in each mouse
            for file_name in file_names[dir_counter]:

                #Load location of recording:
                depth = np.loadtxt(file_dir+file_name+'/depth.txt', dtype=str)
                print file_name, " ", depth

                #Load anesthetic state:
                state = np.loadtxt(file_dir+file_name+'/state.txt', dtype=str)
                if state != state_match: continue

                #Load units from .csv file name; NB: This may lead to duplication as resaved .csv's for Dongsheng
                files = os.listdir(file_dir+file_name)
                temp_names = []
                for file_ in files:
                    if ("unit_" in file_) and ('.csv' in file_):
                        temp_names.append(file_)

                #Save individual unit names, channels, ptps
                units = []
                channels = []
                ptps = []
                for file_ in temp_names:
                    units.append(int(file_[5:7]))
                    channels.append(int(file_[16:18]))
                    ptps.append(int(file_[23:26]))

                for i in range(len(units)):
                    #print "Unit: ", units[i]

                    #Load .csv file to get # spikes
                    csv_file = glob.glob(file_dir + file_name + "/unit_"+str(units[i]).zfill(2) + "*.csv")[0]
                    spikes = np.loadtxt(csv_file, dtype=np.float32)
                    
                    #Load .npy ROI file to get max/min DF/F
                    #Load left side only
                    temp = glob.glob(file_dir + file_name + "/"+file_name+"_matrix_unit_"+str(units[i]).zfill(2) + "_left.npy")
                    if len(temp)==0: continue
                    roi_file_left = glob.glob(file_dir + file_name + "/"+file_name+"_matrix_unit_"+str(units[i]).zfill(2) + "_left.npy")[0]
                    roi_left = np.load(roi_file_left)

                    #Load right side only
                    #roi_file_left = glob.glob(file_dir + file_name + "/"+file_name+"_matrix_unit_"+str(units[i]).zfill(2) + "_right.npy")[0]
                    #roi_left = np.load(roi_file_left)
                    
                    #Plot ROIs
                    if False: plot_ROI(roi_left, roi_right,units[i], area_names, file_name)

                    #Quality control:
                    if True:
                        v_abs1 = max(np.max(roi_left), - np.min(roi_left))
                        #v_abs2 = max(np.max(roi_right), - np.min(roi_right))
                        #global_max = max(v_abs1,v_abs2)
                        if (ptps[i] < 35.0) or (len(spikes)<128) or (v_abs1<0.8): continue
                    
                    else: 
                        if (ptps[i] < 35.0) or (len(spikes)<128): continue
                    
                    #Identify type of connection for max relationship: excitatory vs. inhibitory
                    if (np.max(roi_left)> -np.min(roi_left)):
                        connect_type = 1 #Excitatory connection
                        roi_index = np.where(roi_left==roi_left.max())[0][0] #ROI location of activation
                        #left_tag_cortex_excitatory.append(max_index)
                    else:
                        connect_type = 0 #Inhibitory connection
                        roi_index = np.where(roi_left==roi_left.min())[0][0] #ROI location of activation
                        #left_tag_cortex_inhibitory.append(min_index)
                    
                    if depth == 'cortex':
                        #if channels[i] >10: continue #Skip cortical cells recorded on channel 10 or lower...
                        unrolled = roi_left.ravel()
                        max_unrolled = max(max(unrolled),-min(unrolled))
                        unrolled = unrolled/max_unrolled

                        #Interpolate 
                        x = np.linspace(0, len(unrolled), len(unrolled))
                        y = unrolled
                        f = interpolate.interp1d(x,y)
                        x2 = np.linspace(0, len(unrolled), len(area_names)*200)
                        unrolled = f(x2)
                        total_cortex_roi_left[roi_index][connect_type].append(unrolled)
                        total_cortex_roi_left_cumulative[connect_type].append(unrolled)
                        total_cortex_all.append(unrolled)
                        total_cortex_all_files.append(roi_file_left)
                        total_cortex_all_ch.append(channels[i])

                        cortex_unit_counter+=1 

                    else: 

                        max_index = np.where(roi_left==roi_left.max())[0][0]
                        min_index = np.where(roi_left==roi_left.min())[0][0]
                        unrolled = roi_left.ravel()
                        max_unrolled = max(max(unrolled),-min(unrolled))
                        unrolled = unrolled/max_unrolled
                        
                        #Interpolate 
                        x = np.linspace(0,len(unrolled), len(unrolled))
                        y = unrolled
                        f = interpolate.interp1d(x,y)
                        x2 = np.linspace(0, len(unrolled), len(area_names)*200)
                        unrolled = f(x2)
                        total_subcortical_roi_left[roi_index][connect_type].append(unrolled)
                        
                        total_subcortical_roi_left_cumulative[connect_type].append(unrolled)
                        total_subcortical_all.append(unrolled)
                        total_subcortical_all_files.append(roi_file_left)
                        total_subcortical_all_ch.append(channels[i])

                        subcortical_unit_counter+=1 

            mouse_cortex_counter.append(cortex_unit_counter)
            mouse_subcortical_counter.append(subcortical_unit_counter)
            
            dir_counter+=1

        data.cortex = np.array(total_cortex_all)
        data.cortex_files = total_cortex_all_files
        data.cortex_counter = np.array(mouse_cortex_counter)
        data.cortex_ch = total_cortex_all_ch

        np.save(main_dir+"total_cortex_all_"+state_match, data.cortex)
        np.save(main_dir+"total_cortex_all_files_"+state_match, data.cortex_files)
        np.save(main_dir+"total_cortex_all_counter_"+state_match, data.cortex_counter)
        np.save(main_dir+"total_cortex_all_ch_"+state_match, data.cortex_ch)
        
        
        data.subcortical = np.array(total_subcortical_all)
        data.subcortical_files = total_subcortical_all_files
        data.subcortical_counter = np.array(mouse_subcortical_counter)
        data.subcortical_ch = total_subcortical_all_ch
        
        np.save(main_dir+"total_subcortical_all_"+state_match, data.subcortical)
        np.save(main_dir+"total_subcortical_all_files_"+state_match, data.subcortical_files)
        np.save(main_dir+"total_subcortical_all_counter_"+state_match, data.subcortical_counter)
        np.save(main_dir+"total_subcortical_all_ch_"+state_match, data.subcortical_ch)

    else:
        
        data.cortex = np.load(main_dir + "total_cortex_all_"+state_match+'.npy')
        data.cortex_files = np.load(main_dir+"total_cortex_all_files_"+state_match+'.npy')
        data.cortex_counter = np.load(main_dir+"total_cortex_all_counter_"+state_match+'.npy')
        data.cortex_ch = np.load(main_dir+"total_cortex_all_ch_"+state_match+'.npy')

        data.subcortical = np.load(main_dir+"total_subcortical_all_"+state_match+'.npy')
        data.subcortical_files = np.load(main_dir+"total_subcortical_all_files_"+state_match+'.npy')
        data.subcortical_counter = np.load(main_dir+"total_subcortical_all_counter_"+state_match+'.npy')
        data.subcortical_ch = np.load(main_dir+"total_subcortical_all_ch_"+state_match+'.npy')


def plot_clusters(stmtd, n_areas, window, area_names, maxmin, hemisphere, file_dirs, file_names, main_dir):


    n_interp_points = 600
   
    colors=['blue','green','violet','lightseagreen','lightsalmon','violet','dodgerblue','mediumvioletred','indianred','lightsalmon','pink','darkolivegreen']
    subcolors=[]
    for i in range(4):
        subcolors.append(['blue','green','violet','lightseagreen','lightsalmon','violet','dodgerblue','mediumvioletred','indianred','lightsalmon','pink','darkolivegreen'])
    alt_colors = ['brown','violet','dodgerblue','mediumvioletred','indianred','lightseagreen','lightsalmon','pink','darkolivegreen','black']
    cluster_names = ['cortex_awake', 'cortex_anesth', 'subcortical_awake', 'subcortical_anesth']
    area_loop = np.arange(len(area_names))

    #***************************** TRACE CLASSIFICATION AND DISPLAY *************

    saved_membership=[]
    saved_peaks=[]
    saved_troughs=[]
    saved_peaks_indexes=[]
    saved_troughs_indexes=[]
    
    cell_type_array = []
    cell_type_scatter = []

    for kk in [1]:
        saved_templates = []
        
        cluster_traces_saved = []
        xx = np.linspace(-float(window),float(window),n_interp_points)

        #************ LOAD EXISTING CLUSTER CLASSES **************
        if (os.path.exists(main_dir+cluster_names[kk]+"_class_0.npy")==True):
            final_classes= []
            #csv_files = glob.glob(main_dir+cluster_names[kk]+"*.npy")
            csv_files = glob.glob(main_dir+"subcortical_anesth_class_"+"*.npy")

            for j in range(len(csv_files)):
                #final_classes.append(np.load(main_dir+cluster_names[kk]+"_class_"+str(j)+'.npy'))
                final_classes.append(np.load(main_dir+"subcortical_anesth_class_"+str(j)+'.npy'))
            
            n_classes = len(csv_files)
            print "Number of loaded classes: ", n_classes
            
        else:
            pass

        #********************************************************************************************************
        #********************************** LOAD EXISTING CLASSES ************************************************
        print "Refitting traces to final classes..."

        locations = [area_names[2]]  #areas: [ 'hindlimb', 'forelimb', 'barrel', 'motor', 'visual', 'retrosplenial', 'acc', 'allcortex']
        hemispheres = ["left"]
        maxormins = ['max']

        cluster_traces_saved = []  #REset clustering lists
        cluster_traces_saved_nonnorm = []
        channels = []
        rasters = []
        templates=[]
        mouse_id = []
        for s in range(len(locations)):
            cluster_traces_saved.append([])
            cluster_traces_saved_nonnorm.append([])
            channels.append([])
            rasters.append([])
            templates.append([])
            mouse_id.append([])
        
        counter=0
        saved_peaks=[]
        saved_troughs=[]
        saved_peaks_indexes=[]
        saved_troughs_indexes=[]
        #Loop over conditions - for particular state/depth (e.g. cortex-awake);
        for location, hemisphere, maxormin in zip(locations,hemispheres,maxormins):
            aa = location  #Use aa and kk legacy from older code
            state =kk

            cluster_traces_saved[counter], channels[counter], rasters[counter], templates[counter], mouse_id[counter], file_indexes_out = \
            Classify_traces(location, hemisphere, maxormin, state, stmtd.tc_data, stmtd.area_maxormin, stmtd.area_hemisphere, stmtd.depth_left, stmtd.state_left, stmtd.channel_left, 
            stmtd.raster_left, stmtd.template_left, stmtd.area_index_left, stmtd.mouse_index, stmtd.file_index)
            
            #print "number of file_names ", len(file_indexes_out)
            
            cluster_traces_saved_nonnorm[counter]=np.array(cluster_traces_saved[counter]).copy()
            
            #Save clustered trace data for offline PCA analysis:
            #Normalize first
            for trace in range(len(cluster_traces_saved[counter])):
                temp1 = np.array(cluster_traces_saved[counter][trace])   #Take each cell trace
                temp1 = temp1/max(np.max(temp1),abs(np.min(temp1))) #and normalize (divide by largest of max or min values)
                cluster_traces_saved[counter][trace]=temp1
                
            np.savetxt('/media/cat/12TB/in_vivo/tim/dongsheng/clusters/state_'+str(kk),cluster_traces_saved[counter])
            #quit()
                
            #quit()

            #************* Match each trace to loaded classes
            clusters = []  #REset cluster lists
            clusters_nonnorm = []
            cell_membership = [] #Keep track of what cluster/class each cell goes to
            names_files=[]
            for s in range(n_classes):
                clusters.append([])
                clusters_nonnorm.append([])
                cell_membership.append([])
                names_files.append([])

            for trace in range(len(cluster_traces_saved[counter])):
                #Compute the closest matching template - no jitter initially
                temp=[]
                #temp_maxdiff=[]
                for c in range(n_classes): #Loop over existing classes
                    if (kk<2) and (c>2): continue #Force cortex shapes only into first 3 classes; required in order to share classification
                    if len(final_classes[c])>0:
                        #Normalize traces before fit
                        temp1 = np.array(cluster_traces_saved[counter][trace])   #Take each cell trace
                        temp1 = temp1/max(np.max(temp1),abs(np.min(temp1))) #and normalize (divide by largest of max or min values)
                        cluster_traces_saved[counter][trace]=temp1
                        
                        #Normalize classes
                        temp2 = np.array(final_classes[c])
                        temp2 = temp2/max(max(temp2),abs(min(temp2))) #Divide by largest of max or min values
                        #temp_maxdiff.append(max(temp1-temp2))
                        
                        #Jitter value +/- 2 time steps
                        jitter_values = []
                        for jj in range(-2,2,1):
                            temp3 = np.roll(temp1,jj)
                            jitter_values.append(np.sum(np.sqrt(np.square(temp3-temp2))))
                        temp.append(min(np.array(jitter_values)))  #Assign the membership based on the lwest rms fit among the jittered values
                        
                if (min(np.array(temp))<150):# and (max(temp_maxdiff)<1.8):
                    clusters[np.argmin(temp)].append(cluster_traces_saved[counter][trace])
                    cell_membership[np.argmin(temp)].append(trace)

                    clusters_nonnorm[np.argmin(temp)].append(cluster_traces_saved_nonnorm[counter][trace])  #Save for DF/F histograms;
                    names_files[np.argmin(temp)].append(file_indexes_out[trace])
            
            
            ##Save clustered trace data for offline PCA analysis:
            #for c in range(n_classes):
            #    np.savetxt('/media/cat/12TB/in_vivo/tim/dongsheng/clusters/cluster_state_'+str(kk)+"_class_"+str(c),clusters[c])
            #    np.savetxt('/media/cat/12TB/in_vivo/tim/dongsheng/clusters/cluster_state_'+str(kk)+"_class_"+str(c)+"_index",cell_membership[c])
            #continue
            
            #Remove "bleaching-like" artifacts from subcortical recs; not 100% convinced they are all artifacts
            n_bleached_cells=0
            #Remove last cluster from each depth/state classification as it's used for gathering bleaching artifacts
            n_bleached_cells = len(clusters[-1])
            del cell_membership[-1] 
            del clusters[-1]
            del clusters_nonnorm[-1]

            for n in range(5):
                #for i in range(438):
                #    if cluster_labels[i]==n: #If cluster label value matches current class 
                #        if location[i]=='cortex':
                #            class_pca_cortex.append(clusters[i])
                #            class_pca_cortex_index.append(cell_indexes[i])  #Save original index of trace 
                #        else:
                #            class_pca_subcortical.append(clusters[i])
                #            class_pca_subcortical_index.append(cell_indexes[i])
                
                clusters[n] = np.loadtxt('/media/cat/12TB/in_vivo/tim/dongsheng/clusters/classout_'+str(n)+'_state_'+str(kk))
                cell_membership[n] = np.int32(np.loadtxt('/media/cat/12TB/in_vivo/tim/dongsheng/clusters/classout_'+str(n)+'_state_'+str(kk)+'_index'))

            
            print "File names for cluster 0"
            for n in range(5):
                print ""
                print ""
                print "Cluster: ", n
                for m in range(len(cell_membership[n])):
                    #print cell_membership[n][m]
                    print file_indexes_out[cell_membership[n][m]]


                #np.savetxt('/media/cat/12TB/in_vivo/tim/clusters/classout_'+str(n)+'_state_2',class_pca_subcortical)
                #np.savetxt('/media/cat/12TB/in_vivo/tim/clusters/classout_'+str(n)+'_state_2_index',class_pca_subcortical_index)


            #************************** CLUSTER LOOP *******************
            n_plots=9
            n_cols = n_plots
            n_rows=5
            font_size=15
            fwhm_data = []
            if len(cluster_traces_saved[counter])==0: continue

            print "Plotting ..."
            print cell_membership
            saved_membership.append(cell_membership) #; continue
            
            fwhm_x=[]
            fwhm_y=[]
            scatter_color=[]
            for p in range(len(cell_membership)):

                #**********Plot time course clusters + std
                ax = plt.subplot(n_rows,n_plots,counter*n_plots+p+1)
                temp_ave = []
                xx = np.linspace(-float(window),float(window),n_interp_points)
                for pp in range(len(clusters[p])):
                    #print p, pp
                    if len(clusters[p][pp])>0:
                        temp_ave.append(clusters[p][pp])
                        #if p==4:
                        #plt.close()
                        #if clusters[p][pp][300]>0:
                        #    plt.plot(xx,clusters[p][pp], color='black', alpha=0.2)
                        #plt.show()
                        #if (p==4) and (pp == 0): 
                        #    np.save("/media/cat/12TB/in_vivo/tim/new_class.npy",clusters[p][pp])
                        #    #quit()

                xx = np.linspace(-float(window),float(window),n_interp_points)            
                if len(temp_ave)>0:
                    temp_std = np.std(np.array(temp_ave), axis=0)
                    temp_ave = np.average(np.array(temp_ave),axis=0)
                    plt.plot(xx, temp_ave, color=subcolors[kk][p], linewidth=4)
                    ax.fill_between(xx, temp_ave+temp_std, temp_ave-temp_std, facecolor=subcolors[kk][p], alpha=0.3)
                
                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)
                if p == 0:
                    ax.get_yaxis().set_visible(True)
                    n_cells = len(list(itertools.chain.from_iterable(clusters)))
                    plt.ylabel(location+" "+hemisphere+"\n #cells:"+str(n_cells)+"  " + str(int(float(n_cells)/(len(cluster_traces_saved[0])-n_bleached_cells)*100.))+"%",fontsize=font_size)
                    ax.tick_params(axis='y', which='major', labelsize=20)
                    ax.yaxis.set_ticks([])
                plt.title("#cells: "+str(len(clusters[p])), fontsize=12)
                plt.plot([-window,window], [0,0], color='black')
                plt.plot([0,0], [-3,7], color='black')
                plt.ylim(-2,2) 
                plt.xlim(-3,3)

                #**********************Shade excitatory and inhibitory periods in graphs
                if False:
                    if len(temp_ave)>0:
                        exc_period=[]
                        inh_period=[]
                        
                        start_time= 200
                        thresholds =0.05
                        updown='baseline'
                        if temp_ave[start_time]<-thresholds:
                            updown='up'
                            exc_period.append(start_time)
                        elif temp_ave[start_time]>thresholds:
                            updown='down'
                            inh_period.append(start_time)
                        
                        for i in range(start_time,400,5):
                            #check transition to inhibition
                            #if p==2: print i, updown
                            if ((updown=='up') or (updown=='baseline')) and (temp_ave[i]<-thresholds):
                                if updown=='up': exc_period.append(i)
                                updown='down'
                                inh_period.append(i)
                            elif ((updown=='down') or (updown=='baseline')) and (temp_ave[i]>thresholds):
                                if updown=='down': inh_period.append(i)
                                updown='up'
                                exc_period.append(i)
                            elif (temp_ave[i]<thresholds) and (temp_ave[i]>-thresholds):
                                if updown=='up': exc_period.append(i)
                                elif updown=='down': inh_period.append(i)
                                updown='baseline'
                            
                        #Add end points in case odd number of transitions
                        inh_period.append(400)
                        inh_period = np.array(inh_period)/100. - 3.
                        exc_period.append(400)
                        exc_period = np.array(exc_period)/100. - 3.
                        
                        for i in range(len(inh_period)/2):
                            ax.axvspan(inh_period[i*2], inh_period[i*2+1], facecolor='grey', alpha=0.4)
                        
                        for i in range(len(exc_period)/2):
                            ax.axvspan(exc_period[i*2], exc_period[i*2+1], facecolor='red', alpha=0.2)
                    
                    #ax.fill_between(xx, temp_ave+temp_std, temp_ave-temp_std, facecolor=subcolors[kk][p], alpha=0.3)
                
                
                #**********Compute template averages and FWHM cell typing structures
                if False:
                    temp_templates = np.float32(templates[counter])[cell_membership[p]]
                    temp_ave2 =[]
                    n_interp_points=600
                    if len(cell_membership[p])>0:
                        for cm in range(len(cell_membership[p])):
                            x_int = np.arange(0, len(temp_templates[cm]), 1)
                            f = interpolate.interp1d(x_int, temp_templates[cm])
                            xx = np.linspace(x_int[0],x_int[-1],n_interp_points)
                            data_ = f(xx)/max(np.abs(f(xx)))
                            #Search for fwhm of trough
                            for t in range(len(xx)):
                                if float(data_[t])<-0.5: x1=t; break
                            for t in range(len(xx)-1,0,-1):
                                if float(data_[t])<-0.5: x2=t; break
                            y1 = float(x2-x1)/7.5/20. #Convert to ms
                            #Search for fwhm of peak
                            for t in range(300,len(xx),1): #Search from midpoint to end
                                if float(data_[t])>(np.max(data_[300:])*.5): x1=t; break
                            #x1 = np.argmin(data_)
                            for t in range(len(xx)-1,300,-1): #Search from midpoint to end
                                if float(data_[t])>(np.max(data_[300:])*.5): x2=t; break
                            #x2 = np.argmax(data_[300:])
                            y2 = float(x2-x1)/7.5/20.
                            #print y1,y2, len(fwhm_x)
                            fwhm_x.append(y1)
                            fwhm_y.append(y2)
                            scatter_color.append(p) #preserve color 
                            cell_type_scatter.append([y1,y2])

                            #Compute average shape + std
                            temp_plot = np.roll(temp_templates[cm],40-np.argmin(temp_templates[cm]))/max(np.abs(temp_templates[cm])) #Line up templates to min value
                            cell_type_array.append(temp_plot)
                            temp_ave2.append(temp_plot)
                        
                        
                        temp_ave2 = np.array(temp_ave2)[:,30:70]
                        temp_std = np.std(np.array(temp_ave2), axis=0)
                        temp_ave2 = np.average(np.array(temp_ave2),axis=0)/1.1+1.5
                        saved_spike_ave = temp_ave2
                        xy=np.arange(-10,30,1.)/13.+1.1
                        print len(xy),len(temp_ave2)
                        if True: #Skip plotting of cell shapes
                            plt.plot(xy,temp_ave2, color=subcolors[kk][p], linewidth=1)
                            ax.fill_between(xy, temp_ave2+temp_std, temp_ave2-temp_std, facecolor=subcolors[kk][p], alpha=0.4)

                #**********Plot channel distribution
                if False:
                    ax = plt.subplot(n_rows,n_plots,(counter+2)*n_plots+p+1)
                    temp_channels = np.array(channels[counter])[cell_membership[p]]
                    xy = np.arange(0, 1.01, .1)
                    yy =  np.histogram(temp_channels, bins = xy)
                    plt.barh(yy[1][:-1],yy[0],.1,color=subcolors[kk][p])
                    plt.ylim(1,0)
                
                #**********Plot isi distributions
                if False:
                    ax = plt.subplot(n_rows,n_plots,(counter+3)*n_plots+p+1)
                    temp_isi = []
                    temp_isi2 = []
                    f_rate=[]
                    for cm in range(len(cell_membership[p])):
                        raster_temp = np.array(rasters[counter])[cell_membership[p]][cm] #Select raster from the state/depth separated and class separated lists
                        for q in range(1,len(raster_temp)-2,1):
                            pre_spike = (raster_temp[q]-raster_temp[q-1])/Sort.samplerate*1E3 #Convert to ms time
                            post_spike = (raster_temp[q+1]-raster_temp[q])/Sort.samplerate*1E3
                            if pre_spike == 0: continue     #Duplicates in data
                            if post_spike == 0: continue    #Duplicates in data
                            #if pre_spike < 0. or post_spike <0: print pre_spike, post_spike; quit()
                            temp_isi.append(pre_spike)
                            temp_isi2.append(post_spike)
                        f_rate.append(len(raster_temp)/(raster_temp[-1]/Sort.samplerate))
                    
                    if len(temp_isi)>0: 
                        print len(cell_membership[p])
                        plt.scatter(np.array(temp_isi), np.array(temp_isi2), s=4, alpha=(.526/(len(cell_membership[p])*2)))#, color=colors[labels])
                        f_rate_std = np.std(np.array(f_rate), axis=0)
                        f_rate=np.average(np.array(f_rate),axis=0)
                        plt.title(str(int(f_rate)) + "+/-" +str(int(f_rate_std))+"hz")
                        ax.set_xscale('log')
                        ax.set_yscale('log')
                        ax.set_xlim(1E0,2E3)
                        ax.set_ylim(1E0,2E3)

                #****************Plot mouse distribution pie charts
                if False:
                    ax = plt.subplot(n_rows,n_plots,(counter+1)*n_plots+p+1)
                    mouse_temp = np.array(mouse_id[counter])[cell_membership[p]]
                    aaa = []
                    for p in np.unique(mouse_temp):
                        a_temp = len(np.where(mouse_temp==p)[0])
                        #print p, a_temp
                        aaa.append(a_temp)
                    aaa = np.float32(aaa)
                    aaa = aaa/np.sum(aaa)
                    #print aaa
                    plt.xlim(-0.05,0.05)
                    plt.ylim(-0.05,0.05)
                    ax.get_yaxis().set_visible(False)
                    ax.get_xaxis().set_visible(False)
                    temp_colors = np.array(alt_colors)[np.unique(mouse_temp)]
                    draw_pie(ax,aaa, 0.0, 0.0, size=10000, colors=temp_colors) #Preserve mouse # and color even when a mouse not present in pie chart

                #**************************** PLOT METRICS AT END OF CLASSIFICATION CHARTS ********************
                if False:
                    #counter = 1
                    tick_size = 10
                    font_size=12
                    cortex_peaks = []
                    cortex_troughs = []
                    cortex_peaks_indexes = []
                    cortex_troughs_indexes= []

                    #Cortex data: well behaved, no need to apply heuristics
                    #Find max peak time first
                    peaks=[]
                    troughs=[]
                    peaks_indexes=[]
                    troughs_indexes = []
                    traces_saved = np.array(clusters_nonnorm[p])
                    for t in range(len(traces_saved)):
                        #print np.min(cluster_traces_saved[t])                #DF/F values
                        
                        temp_peak = np.argmax(traces_saved[t,250:])+250       #Time index of peaks; only look for peaks around t = 0 - otherwise likely poor quality traces
                        temp_trough = np.argmin(traces_saved[t,250:])+250       #Time index of troughs: troughs from t=-1 .. 3.0 sec

                        #if (temp >250):# and (temp < 400):
                        if (traces_saved[t][temp_peak])>abs(traces_saved[t][temp_trough]):
                            peaks.append(traces_saved[t][temp_peak])                #DF/F values
                            peaks_indexes.append(temp_peak)     
                        else:
                            troughs.append(abs(traces_saved[t][temp_trough]))
                            troughs_indexes.append(temp_trough)

                    
                    #PLOT PEAK TIME HISTOGRAMS
                    ax = plt.subplot(n_rows,n_plots,(counter+1)*n_plots+p+1)
                    if p <2: 
                        plt.xlim(-1,1)
                        width = 5
                    else:
                        plt.xlim(-3,3)
                        width=15
                    xy = np.arange(0,600,width)
                    yy =  np.histogram(peaks_indexes, bins = xy)
                    xx = (yy[1][:-1]-300)/100.
                    plt.bar(xx,yy[0],width/100., color=colors[p])
                    temp1=max(yy[0])
                    if p == 0: plt.ylabel("Time-to-peaks \n(x-axis: sec)",fontsize=font_size)
                    ax.tick_params(axis='both', which='both', labelsize=tick_size)
                    plt.plot([-3,3],[0,0],color='black',linewidth=1)
                    ax.xaxis.set_ticks(np.arange(-3.0,3.1,1.))
                    ax.yaxis.set_ticks([])
                    cortex_peaks_indexes.append(peaks_indexes)
                    plt.plot([0,(i+1)*3],[0,0],color='black',linewidth=1)

                    #PLOT TROUGH TIME HISTOGRAMS
                    xy = np.arange(0,600,width)
                    yy =  np.histogram(troughs_indexes, bins = xy)
                    xx = (yy[1][:-1]-300)/100.
                    if ((kk!=0) or (p!=0)) and ((kk!=2) or (p!=1)):
                        plt.bar(xx,-yy[0],width/100., color=colors[p],hatch = '//')    
                        cortex_troughs_indexes.append(troughs_indexes)

                    temp2=max(yy[0])
                    ax.tick_params(axis='both', which='both', labelsize=tick_size)
                    ax.yaxis.set_ticks([])
                    ax.xaxis.set_ticks(np.arange(-3.0,3.1,1.))
                    max_lim = max(temp1,temp2)
                    plt.ylim(-max_lim,max_lim)
                    
                    #PLOT PEAK DF/F
                    ax = plt.subplot(n_rows,n_plots,(counter+2)*n_plots+p+1)
                    width = .25
                    xy = np.arange(0,10,.25 )
                    yy =  np.histogram(peaks, bins = xy)
                    plt.bar(yy[1][:-1],yy[0],width, color=colors[p])
                    cortex_peaks.append(peaks)
                    temp1=max(yy[0])
                    if counter%n_cols ==1:  plt.ylabel(location+ " " + hemisphere + "\n#cells: " +str(len(peaks)))
                    if p == 0: plt.ylabel("DF/F - Peaks\n(x-axis: %DF/F)",fontsize=font_size)
                    ax.yaxis.set_ticks([])
                    plt.plot([-3,3],[0,0],color='black',linewidth=1)
                    ax.tick_params(axis='both', which='both', labelsize=tick_size)
                    plt.xlim(0,6)
                    
                    #PLOT TROUGH DF/F
                    width = .25
                    xy = np.arange(0,10,.25 )
                    yy =  np.histogram(troughs, bins = xy)
                    if ((kk!=0) or (p!=0)) and ((kk!=2) or (p!=1)):
                        plt.bar(yy[1][:-1],-yy[0],width, color=colors[p], hatch = '//')
                        cortex_troughs.append(troughs)
                    
                    temp2=max(yy[0])
                    ax.yaxis.set_ticks([])

                    ax.tick_params(axis='both', which='both', labelsize=tick_size)
                    plt.xlim(0,6)

                    max_lim = max(temp1,temp2)
                    plt.ylim(-max_lim,max_lim)
                    
                    #CONVERT LISTS OF PEAKS/TROUGHS and DF/F VALUES TO NUMPY ARRAYS
                    saved_peaks.append(np.array(cortex_peaks))
                    saved_troughs.append(np.array(cortex_troughs))
                    saved_peaks_indexes.append((np.array(cortex_peaks_indexes)-300)/100.)
                    saved_troughs_indexes.append((np.array(cortex_troughs_indexes)-300)/100.)
                    
            #**********Plot scatter FWHM distributions
            ax = plt.subplot(n_rows,n_plots,(counter)*n_plots+n_plots-2)

            if False:
                for i in range(len(fwhm_x)):
                    plt.scatter(fwhm_x[i],fwhm_y[i], s=2, color=subcolors[0][scatter_color[i]]) #subcolors[kk][pp])
                xx=np.linspace(0.45,0.75,len(saved_spike_ave))
                plt.plot(xx, saved_spike_ave/1.5-.3, color='black', linewidth=2)
                plt.xlim(0,.8)
                ax.xaxis.set_ticks(np.arange(0,.81,.2))
                plt.ylim(0,1.5)
                plt.tick_params(axis='both', which='major', labelsize=10)
                plt.title("Putative Cell Clasification", fontsize=font_size)
                plt.xlabel("Trough FWHM (ms)", fontsize=font_size)
                plt.ylabel("Peak FWHM (ms)", fontsize=font_size)

            #Plot class distribution pie charts
            ax = plt.subplot(n_rows,n_plots,counter*n_plots+n_plots-1)
            a = []
            for p in range(len(cell_membership)):
                a.append(len(cell_membership[p]))
            a = np.float32(a)
            a = a/np.sum(a)

            draw_pie(ax,a, 0.0, 0.0, size=6000, colors=subcolors[kk])
            plt.xlim(-0.05,0.05)
            plt.ylim(-0.05,0.05)
            ax.get_yaxis().set_visible(False)
            ax.get_xaxis().set_visible(False)

            counter+=1 

            #**********************Distributions plots of DF/F and time-to-peak
            if kk==0: 
                temp_clusters = 3
            else: 
                temp_clusters=5
            if False:
                #Plot time-to-peak distributions
                ax = plt.subplot(n_rows,n_plots,n_plots+6)
                maxtemp=0
                xx = np.arange(0,99,1)
                for i in range(temp_clusters):
                    temp1 = np.average(saved_peaks_indexes[i])
                    temp2 = np.average(saved_troughs_indexes[i])
                    ax.barh(xx[i*2],temp1,1,color=colors[i],alpha=1.0)
                    ax.barh(xx[i*2]+1,temp2,1,color=colors[i], hatch = '//', alpha=1.0)
                    ax.errorbar(temp1, xx[i*2]+.5,xerr = scipy.stats.sem(saved_peaks_indexes[i][0]), color='black', linewidth = 3, capsize=6, capthick=1, alpha=1)
                    if ((kk!=0) or (i!=0)) and ((kk!=2) or (i!=1)): ax.errorbar(temp2, xx[i*2]+1.5,xerr = scipy.stats.sem(saved_troughs_indexes[i][0]), color='black', linewidth = 3, capsize=6, capthick=1, alpha=1)
                ax.get_yaxis().set_visible(False)
                plt.ylim(xx[0],xx[i*2]+3)
                plt.xlim(-3,3)
                plt.plot([0,0],[xx[0],xx[-1]],color='black',linewidth=1)
                ax.tick_params(axis='both', which='both', labelsize=tick_size)

                #Plot DF/F peaks distributions
                xx = [1,4,7,10,13,16,19,22]
                ax = plt.subplot(n_rows,n_plots,n_plots*2+6)
                maxtemp=0
                for i in range(temp_clusters):
                    temp1 = np.average(saved_peaks[i])
                    temp2 = -np.average(saved_troughs[i])
                    ax.bar(xx[i],temp1,1.5,color=colors[i],alpha=1.0)
                    ax.bar(xx[i],temp2,1.5,color=colors[i], hatch = '//', alpha=1.0)
                    ax.errorbar(xx[i]+.75, temp1,yerr = scipy.stats.sem(saved_peaks[i][0]), color='black', linewidth = 6, capsize=12, capthick=2, alpha=1)
                    if ((kk!=0) or (i!=0)) and ((kk!=2) or (i!=1)): ax.errorbar(xx[i]+.75, temp2,yerr = scipy.stats.sem(saved_troughs[i][0]), color='black', linewidth = 6, capsize=10, capthick=2, alpha=1)
                ax.get_xaxis().set_visible(False)
                plt.xlim(0,(i+1)*3)
                plt.ylim(-3,3)
                plt.plot([0,(i+1)*3],[0,0],color='black',linewidth=1)
                ax.tick_params(axis='both', which='both', labelsize=tick_size)

        plt.suptitle("Area: "+location+ " "+hemisphere, fontsize = 25)
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()    
        plt.show()

    np.save('/media/cat/12TB/in_vivo/tim/cell_type_array', cell_type_array)
    np.save('/media/cat/12TB/in_vivo/tim/cell_type_scatter' , cell_type_scatter)

def cluster_traces(traces, stmtd, state, depth, n_clusters, marker):

    colors=['blue','green','violet','lightseagreen','lightsalmon','dodgerblue','mediumvioletred','indianred','lightsalmon','pink','darkolivegreen']


    #PCA
    X = traces
    n_components = 3
    pca_data, coords = PCA(X, n_components)

    
    #********* KMEANS CLUSTERING *******
    if True:
        #n_clusters = 5
        cluster_labels = KMEANS(pca_data, n_clusters)
        
    else:
        cluster_labels = Meanshift(pca_data)
        n_clusters = np.max(cluster_labels)+1

    cell_indexes=[]
    color_index = []
    for i in range(len(X)):
        color_index.append(colors[cluster_labels[i]])
        cell_indexes.append(cluster_labels[i])

    #PLOT PCA CLUSTERS
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_title(state+"  " + depth + " #clusters: " + str(n_clusters))
    cmhot = plt.get_cmap("hot")
    
    for i in range(len(cluster_labels)):
        ax.scatter(coords[0][i],coords[1][i],coords[2][i], s=120, c=color_index[i], edgecolor=color_index,marker=marker[i])

    plt.show()

    for k in range(n_clusters):
        ax1=plt.subplot(2,3,k+1)
        plt.plot([0,0],[-1,1], color='black')
        plt.plot([-3,3],[0,0], color='black')
        plt.ylim(-1,1)
        plt.xlim(-3,3)        
        
        ave_trace = []
        for j in range(len(traces)):
            if cluster_labels[j]==k:
                #plt.plot(traces[j], color=colors[k], linewidth=2, alpha=0.3)
                ave_trace.append(traces[j])

        xx = np.linspace(-3.,3.,600)
        temp_std = np.std(np.array(ave_trace), axis=0)
        temp_ave = np.average(np.array(ave_trace),axis=0)
        ax1.fill_between(xx, temp_ave+temp_std, temp_ave-temp_std, facecolor=colors[k], alpha=0.3)
        
        ax1.plot(xx, temp_ave,color=colors[k], linewidth=5, alpha=1.0)
        plt.title("Cluster: "+str(k+1)+ "  #traces: "+str(len(ave_trace)))
    
    plt.show()



def PCA(X, n_components):
    from sklearn import decomposition

    pca = decomposition.PCA(n_components)
    pca.fit(X)
    X=pca.transform(X)

    coords = []
    for i in range(len(X)):
         coords.append([X[i][0], X[i][1], X[i][2]])
    
    return X, np.array(coords).T













