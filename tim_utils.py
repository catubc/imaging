from tsf_ptcs_classes import *
from distributions import *
from sequential_firing import *

from scipy import stats, signal

import numpy as np
import time, math
import sys
import os.path
import multiprocessing as mp
import matplotlib.animation as anim
from matplotlib import animation
from matplotlib.path import Path

from pylab import *
from scipy.interpolate import interp1d
import struct, array, csv
import scipy.optimize 
from scipy.optimize import curve_fit
import pandas as pd
from scipy.interpolate import UnivariateSpline
import matplotlib.mlab as mlab
from PIL import Image

def Load_images(file_dir, file_name):
    
    if (os.path.exists(file_dir + file_name + '/' + file_name + '_images.npy')==False):

        img = Image.open(file_dir + file_name + '/' + file_name+'.tif')

        counter=0
        while True:
            try:
                img.seek(counter)
            except EOFError:
                break
            counter+=1

        #Default pic sizes
        n_pixels = 256

        #Initialize 3D image array
        n_frames = counter
        images_raw = np.zeros((n_frames, n_pixels, n_pixels), dtype = np.float16)

        print "n_frames: ", n_frames
        for i in range(n_frames): 
            try:
                img.seek(i)
                print "Loading frame: ", i
                images_raw [i] = np.float16(img)
                if i>9640:# and i%10==0: # 
                    im = plt.imshow(images_raw[i])
                    plt.title("Frame: " +str(i))
                    plt.show()

            except EOFError:
                break

        images_start= 0
        images_end = 9642
        images_raw=images_raw[images_start:images_end]

        print "Re-saving imaging array..."

        np.save(file_dir + file_name + '/' + file_name + '_images', images_raw)
        np.savetxt(file_dir + file_name + '/' + file_name + '_images_start_'+str(images_start)+'_end_'+
        str(images_end), [images_start, images_end])

        quit()

    else:
        images_raw = np.load(file_dir + file_name + '/' + file_name + '_images.npy')
        images_raw = np.float16(images_raw)
        return images_raw
        

def Load_images_start_end(file_dir, file_name, images_raw):
    
    data = MCD_read_imagingtimes_old(file_dir + file_name+ '/' + file_name + '.mcd')

    temp_data = []
    for i in range(data['extra'].item_count):
        temp_data.append(data['extra'].get_data(i))

    img_start = temp_data[0][0]
    img_end = temp_data[len(temp_data)-1][0]
            
    n_pixels = len(images_raw[0])
    n_frames = len(images_raw)
    len_frame = float(img_end - img_start) / n_frames
    img_rate = 1. / len_frame

    img_times = []
    for i in range(n_frames):
        img_times.append(img_start+i*len_frame)
    img_times = np.array(img_times)
    
    return img_start, img_end, len_frame, img_rate, n_pixels, n_frames, img_times

def Spike_averages((temp3)):
    global images_temp
    return images_temp[temp3]

def Sum_list((temp_list)):
    global temp_window, temp_img_rate, temp_n_pixels
    temp_sum = np.zeros((int(temp_window*temp_img_rate)*2, temp_n_pixels, temp_n_pixels), dtype=np.float16)
    for i in range(len(temp_list)):
        temp_sum += temp_list[i]
    return temp_sum
    
def Compute_spike_triggered_average(unit, channel, spikes, window, img_rate, img_times, n_pixels, images_raw, file_dir, file_name):
    global images_temp, temp_window, temp_img_rate, temp_n_pixels
    temp_window = window
    temp_img_rate = img_rate
    temp_n_pixels = n_pixels

    n_procs = 12
    print "No. of processors: ", n_procs

    #Remove spikes outside of imaging frames
    print "Total no. spikes: ", len(spikes)
    temp0 = np.where(np.logical_and(spikes>=img_times[0]+window, spikes<=img_times[-1]-window))[0]
    spikes = spikes[temp0]
    print "No. spikes within imaging period: ", len(spikes)

    all_spikes = True
    if all_spikes: 
        plot_string = 'all'
    else:
        plot_string = '1sec'

    #Check to see if images already loaded and saved as .npy
    if (os.path.exists(file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'.npy')==False):
        #Remove baseline
        if False:
            images_temp = images_raw.copy()
            baseline = np.mean(images_temp, axis=0, dtype=np.float32)
            images_temp = (images_temp - baseline)/baseline

        #Trigger only off spikes separated by 1.0 secs
        if all_spikes:
            #Use all spikes
            temp3 = []
            for i in spikes:
                temp3.append(np.where(np.logical_and(img_times>=i-window, img_times<=i+window))[0][0:2*int(window*img_rate)])

        else:
            temp4 = []
            temp4.append(spikes[0])
            counter = 0
            for i in range(1,len(spikes),1):
                if (spikes[i]-temp4[counter])>1.0:
                    temp4.append(spikes[i])
                    counter+=1
            print "Spikes w/in 1sec: ", counter
            spikes = temp4

            temp3 = []
            for i in spikes:
                temp3.append(np.where(np.logical_and(img_times>=i-window, img_times<=i+window))[0][0:2*int(window*img_rate)])
                
        #Compute all frames based on image index;
        print "Computing images from spike averages in parallel ..."
        pool = mp.Pool(n_procs)
        images_triggered_temp = pool.map(Spike_averages, temp3)
        pool.close()
        print "... done "        
        
        #Sum over all spikes - MUST FIND FASTER METHOD
        print "Summing frames over all spikes"
        #temp_images = np.sum(images_triggered_temp,axis=0)
        temp_images = np.zeros((int(window*img_rate)*2, n_pixels, n_pixels), dtype=np.float16)
        for i in range(len(images_triggered_temp)):
            temp_images += images_triggered_temp[i]
        print "... done "        

        #Parallelize summing of frames 
        #temp_list = []
        #for i in range(int(len(images_triggered_temp)/10)-1):
            #temp_list.append(images_triggered_temp[i*int(len(images_triggered_temp)/10):(i+1)*int(len(images_triggered_temp)/10)])

        ##Must also append last bit of array in case not divisible by 10
        #if (float(len(images_triggered_temp))/10.0).is_integer():
            #pass
        #else:
            #temp_list.append(images_triggered_temp[int(len(images_triggered_temp)/10)*10:])
        
        #temp66 = np.zeros((int(window*img_rate)*2, n_pixels, n_pixels), dtype=np.float16)

        #print "Computing list sum in parallel ..."
        #pool = mp.Pool(n_procs)
        #images_triggered_temp2 = pool.map(Sum_list, temp_list)
        #pool.close()
        #print "...done"
        
        #temp_images = np.zeros((int(window*img_rate)*2, n_pixels, n_pixels), dtype=np.float16)
        #print len(images_triggered_temp2)
        #for i in range(len(images_triggered_temp2)):
            #temp_images += images_triggered_temp2[i]
        
        print "No. spikes (counter) : ", len(spikes)
        images_processed = temp_images/float(len(spikes))

        np.save(file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2), images_processed)

    else: 
        print "Skipping processing ... loading from file"
        images_processed = np.load(file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+ '.npy')

    #print "Shape mean: ", images_processed.shape
    
    return images_processed, spikes, plot_string

def on_click(event):
    
    global coords, images_temp
    
    if event.inaxes is not None:
        coords.append((event.ydata, event.xdata))
        for j in range(len(coords)):
            for k in range(7):
                for l in range(7):
                    images_temp[100][min(255,int(coords[j][0])-3+k)][min(255,int(coords[j][1])-3+l)]=0

        ax.imshow(images_temp[100])
        ax.set_title("Compute generic (outside the brain) mask")
        plt.show()

    else:
        print 'Exiting'
        plt.close()
        fig.canvas.mpl_disconnect(cid)

def remove_bregma(event):
    
    global bregma_coords, images_temp, n_pix
    
    if event.inaxes is not None:
        bregma_coords.append((event.ydata, event.xdata))
        for j in range(len(bregma_coords)):
            for k in range(n_pix):
                for l in range(7):
                    images_temp[100][k][min(n_pix-1,int(bregma_coords[j][1])-3+l)]=0

        plt.close()
        fig.canvas.mpl_disconnect(cid)

def define_area(event):
    
    global area_coords, images_temp, n_pix
    
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
        print 'Exiting'
        plt.close()
        fig.canvas.mpl_disconnect(cid)
        
def Compute_generic_mask(images_processed, images_raw, file_dir, file_name):

    global coords, images_temp, ax, fig, cid
    images_temp = images_processed.copy()
    fig, ax = plt.subplots()

    coords=[]

    ax.imshow(images_processed[100])
    ax.set_title("Compute generic (outside the brain) mask")
    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    cid = fig.canvas.mpl_connect('button_press_event', on_click)
    plt.show()

    #Search points outside and black them out:
    all_points = []
    for i in range(len(images_raw[0][0])):
        for j in range(len(images_raw[0][0])):
            all_points.append([i,j])

    all_points = np.array(all_points)
    vertixes = np.array(coords) 
    vertixes_path = Path(vertixes)
    
    mask = vertixes_path.contains_points(all_points)
    counter=0
    coords_save=[]
    for i in range(len(images_raw[0][0])):
        for j in range(len(images_raw[0][0])):
            if mask[counter] == False:
                images_processed[100][i][j]=0
                coords_save.append([i,j])
            counter+=1

    fig, ax = plt.subplots()
    ax.imshow(images_processed[100])
    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    plt.show()
   
    genericmask_file = file_dir + 'genericmask.txt'
    np.savetxt(genericmask_file, coords_save)

    print "Finished Making General Mask"


def Define_cortical_areas(unit, channel, images_processed, file_dir, file_name, n_pixels):
    print "Defining cortical areas"

    global coords, bregma_coords, area_coords, images_temp, ax, fig, cid, n_pix

    images_temp = np.array(images_processed).copy()
    n_pix = n_pixels
    
    #Load Generic Mask
    if (os.path.exists(file_dir + 'genericmask.txt')==True):
        print "Loading existing generic mask"
        generic_mask_file = file_dir + 'genericmask.txt'
        coords_generic = np.loadtxt(generic_mask_file)
        #Ablate generic map
        for i in range(len(images_temp)):
            for j in range(len(coords_generic)):
                images_temp[i][min(255,int(coords_generic[j][0]))][min(255,int(coords_generic[j][1]))]=0
                
    #Load existing artifact mask file
    coords=[]
    specific_mask_file = file_dir +'artifactmask.txt'
    if (os.path.exists(specific_mask_file)==True):
        temp_data= np.loadtxt(specific_mask_file)
        for i in range(len(temp_data)):
            coords.append(temp_data[i])
        update_length=len(coords)

        #Ablate specific map
        for i in range(len(images_temp)):
            for j in range(len(coords)):
                for k in range(3):
                    for l in range(3):
                        images_temp[i][min(255,int(coords[j][0])-3+k)][min(255,int(coords[j][1])-3+l)]=0
    else:
        update_length=0

    #******** Define bregma and auto-remove centreline artifacts *********
    bregma_mask_file = file_dir + 'bregmamask.txt'
    if (os.path.exists(bregma_mask_file)==False):
        bregma_coords = []
        fig, ax = plt.subplots()
        ax.imshow(images_temp[100])
        ax.set_title("Define Bregma location")
        cid = fig.canvas.mpl_connect('button_press_event', remove_bregma)
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        plt.show()

        #Save bregma and centreline mask 
        if len(bregma_coords)>0: 
            np.savetxt(bregma_mask_file, np.int16(bregma_coords))
        bregma_coords = bregma_coords[0] #Save data in single tuple
    else:
        bregma_coords = np.loadtxt(bregma_mask_file)
        #Remove centreline artifacts
        for j in range(len(bregma_coords)):
            for k in range(n_pix):
                for l in range(7):
                    images_temp[100][k][min(n_pix-1,int(bregma_coords[1])-3+l)]=0

    #print "bregma_coords: ", bregma_coords
    #************ Define hindlimb - Green ************
    hindlimb_file = file_dir + 'hindlimb.txt'
    if (os.path.exists(hindlimb_file)==False):
        area_coords = []
        fig, ax = plt.subplots()
        ax.imshow(images_temp[100])
        ax.set_title("Define Hindlimb Location")
        cid = fig.canvas.mpl_connect('button_press_event', define_area)
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        plt.show()

        #Convert coords into x and y data for plotting; select each column, then tack on first coordinate to complete circle
        area_coords = np.array(area_coords)
        hindlimb_xx = np.append(area_coords[:,1], area_coords[0][1])
        hindlimb_yy = np.append(area_coords[:,0], area_coords[0][0])

        #Save total map containing specific coords
        if len(area_coords)>0: 
            np.savetxt(hindlimb_file, np.int16(area_coords))
    else:
        hindlimb_coords = np.loadtxt(hindlimb_file)
        hindlimb_xx = np.append(hindlimb_coords[:,1], hindlimb_coords[0][1])
        hindlimb_yy = np.append(hindlimb_coords[:,0], hindlimb_coords[0][0])
    
    #******** Define Barrel Cortex - Yellow ********
    barrel_file = file_dir + 'barrel.txt'
    if (os.path.exists(barrel_file)==False):
        area_coords = []
        fig, ax = plt.subplots()
        ax.imshow(images_temp[100])
        plt.plot(hindlimb_xx,hindlimb_yy,color='green',linewidth=10)
        plt.plot(bregma_coords[1]-(hindlimb_xx-bregma_coords[1]),hindlimb_yy,color='green',linewidth=10)
        plt.xlim(0,n_pixels-1)
        plt.ylim(n_pixels-1,0)
        ax.set_title("Define Barrel Cortex Location")
        cid = fig.canvas.mpl_connect('button_press_event', define_area)
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        plt.show()

        #Convert coords into x and y data for plotting; select each column, then tack on first coordinate to complete circle
        area_coords = np.array(area_coords)
        barrel_yy = np.append(area_coords[:,0], area_coords[0][0])
        barrel_xx = np.append(area_coords[:,1], area_coords[0][1])

        #Save total map containing specific coords
        if len(area_coords)>0: 
            np.savetxt(barrel_file, np.int16(area_coords))
    else:
        barrel_coords = np.loadtxt(barrel_file)
        barrel_xx = np.append(barrel_coords[:,1], barrel_coords[0][1])
        barrel_yy = np.append(barrel_coords[:,0], barrel_coords[0][0])

    #Define Motor - Red
    motor_file = file_dir + 'motor.txt'
    if (os.path.exists(motor_file)==False):
        area_coords = []
        fig, ax = plt.subplots()
        ax.imshow(images_temp[100])
        plt.plot(hindlimb_xx,hindlimb_yy,color='green',linewidth=10)
        plt.plot(bregma_coords[1]-(hindlimb_xx-bregma_coords[1]),hindlimb_yy,color='green',linewidth=10)
        plt.plot(barrel_xx,barrel_yy,color='yellow',linewidth=10)
        plt.plot(bregma_coords[1]-(barrel_xx-bregma_coords[1]),barrel_yy,color='yellow',linewidth=10)        
        plt.xlim(0,n_pixels-1)
        plt.ylim(n_pixels-1,0)
        ax.set_title("Define Motor Cortex Location")
        cid = fig.canvas.mpl_connect('button_press_event', define_area)
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        plt.show()

        #Convert coords into x and y data for plotting; select each column, then tack on first coordinate to complete circle
        area_coords = np.array(area_coords)
        motor_yy = np.append(area_coords[:,0], area_coords[0][0])
        motor_xx = np.append(area_coords[:,1], area_coords[0][1])

        #Save total map containing specific coords
        if len(area_coords)>0: 
            np.savetxt(motor_file, np.int16(area_coords))
    else:
        motor_coords = np.loadtxt(motor_file)
        motor_xx = np.append(motor_coords[:,1], motor_coords[0][1])
        motor_yy = np.append(motor_coords[:,0], motor_coords[0][0])
        
    #Define Visual - Blue
    visual_file = file_dir + 'visual.txt'
    if (os.path.exists(visual_file)==False):
        area_coords = []
        fig, ax = plt.subplots()
        ax.imshow(images_temp[100])
        plt.plot(hindlimb_xx,hindlimb_yy,color='green',linewidth=10)
        plt.plot(bregma_coords[1]-(hindlimb_xx-bregma_coords[1]),hindlimb_yy,color='green',linewidth=10)
        plt.plot(barrel_xx,barrel_yy,color='yellow',linewidth=10)
        plt.plot(bregma_coords[1]-(barrel_xx-bregma_coords[1]),barrel_yy,color='yellow',linewidth=10)      
        plt.plot(motor_xx,motor_yy,color='red',linewidth=10)
        plt.plot(bregma_coords[1]-(motor_xx-bregma_coords[1]),motor_yy,color='red',linewidth=10)
        plt.xlim(0,n_pixels-1)
        plt.ylim(n_pixels-1,0)
        ax.set_title("Define Visual Cortex Location")
        cid = fig.canvas.mpl_connect('button_press_event', define_area)
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        plt.show()

        #Convert coords into x and y data for plotting; select each column, then tack on first coordinate to complete circle
        area_coords = np.array(area_coords)
        visual_yy = np.append(area_coords[:,0], area_coords[0][0])
        visual_xx = np.append(area_coords[:,1], area_coords[0][1])

        #Save total map containing specific coords
        if len(area_coords)>0: 
            np.savetxt(visual_file, np.int16(area_coords))
    else:
        visual_coords = np.loadtxt(visual_file)
        visual_xx = np.append(visual_coords[:,1], visual_coords[0][1])
        visual_yy = np.append(visual_coords[:,0], visual_coords[0][0])

    #Define Retrosplenial - Cyan
    rs_file = file_dir + 'rs.txt'
    if (os.path.exists(rs_file)==False):
        area_coords = []
        fig, ax = plt.subplots()
        ax.imshow(images_temp[100])
        plt.plot(hindlimb_xx,hindlimb_yy,color='green',linewidth=10)
        plt.plot(bregma_coords[1]-(hindlimb_xx-bregma_coords[1]),hindlimb_yy,color='green',linewidth=10)
        plt.plot(barrel_xx,barrel_yy,color='yellow',linewidth=10)
        plt.plot(bregma_coords[1]-(barrel_xx-bregma_coords[1]),barrel_yy,color='yellow',linewidth=10)      
        plt.plot(motor_xx,motor_yy,color='red',linewidth=10)
        plt.plot(bregma_coords[1]-(motor_xx-bregma_coords[1]),motor_yy,color='red',linewidth=10)
        plt.plot(visual_xx,visual_yy,color='cyan',linewidth=10)
        plt.plot(bregma_coords[1]-(visual_xx-bregma_coords[1]),visual_yy,color='cyan',linewidth=10)
        plt.xlim(0,n_pixels-1)
        plt.ylim(n_pixels-1,0)
        ax.set_title("Define Restrosplenial Cortex Location")
        cid = fig.canvas.mpl_connect('button_press_event', define_area)
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        plt.show()

        #Convert coords into x and y data for plotting; select each column, then tack on first coordinate to complete circle
        area_coords = np.array(area_coords)
        rs_yy = np.append(area_coords[:,0], area_coords[0][0])
        rs_xx = np.append(area_coords[:,1], area_coords[0][1])

        #Save total map containing specific coords
        if len(area_coords)>0: 
            np.savetxt(rs_file, np.int16(area_coords))
    else:
        rs_coords = np.loadtxt(rs_file)
        rs_xx = np.append(rs_coords[:,1], rs_coords[0][1])
        rs_yy = np.append(rs_coords[:,0], rs_coords[0][0])    
    
    #Define Anterior Singulate - Magenta
    acc_file = file_dir + 'acc.txt'
    if (os.path.exists(acc_file)==False):
        area_coords = []
        fig, ax = plt.subplots()
        ax.imshow(images_temp[100])
        plt.plot(hindlimb_xx,hindlimb_yy,color='green',linewidth=10)
        plt.plot(bregma_coords[1]-(hindlimb_xx-bregma_coords[1]),hindlimb_yy,color='green',linewidth=10)
        plt.plot(barrel_xx,barrel_yy,color='yellow',linewidth=10)
        plt.plot(bregma_coords[1]-(barrel_xx-bregma_coords[1]),barrel_yy,color='yellow',linewidth=10)      
        plt.plot(motor_xx,motor_yy,color='red',linewidth=10)
        plt.plot(bregma_coords[1]-(motor_xx-bregma_coords[1]),motor_yy,color='red',linewidth=10)
        plt.plot(visual_xx,visual_yy,color='blue',linewidth=10)
        plt.plot(bregma_coords[1]-(visual_xx-bregma_coords[1]),visual_yy,color='blue',linewidth=10)
        plt.plot(rs_xx,rs_yy,color='cyan',linewidth=10)
        plt.plot(bregma_coords[1]-(rs_xx-bregma_coords[1]),rs_yy,color='cyan',linewidth=10)
        
        plt.xlim(0,n_pixels-1)
        plt.ylim(n_pixels-1,0)
        ax.set_title("Define Anterior Cingulate Cortex Location")
        cid = fig.canvas.mpl_connect('button_press_event', define_area)
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        plt.show()

        #Convert coords into x and y data for plotting; select each column, then tack on first coordinate to complete circle
        area_coords = np.array(area_coords)
        acc_yy = np.append(area_coords[:,0], area_coords[0][0])
        acc_xx = np.append(area_coords[:,1], area_coords[0][1])

        #Save total map containing specific coords
        if len(area_coords)>0: 
            np.savetxt(acc_file, np.int16(area_coords))
    else:
        acc_coords = np.loadtxt(acc_file)
        acc_xx = np.append(acc_coords[:,1], acc_coords[0][1])
        acc_yy = np.append(acc_coords[:,0], acc_coords[0][0])    


    #Plot all parcellated areas
    if False:
        fig, ax = plt.subplots()
        ax.imshow(images_temp[100])
        plt.plot(hindlimb_xx,hindlimb_yy,color='green',linewidth=10)
        plt.plot(bregma_coords[1]-(hindlimb_xx-bregma_coords[1]),hindlimb_yy,color='green',linewidth=10)
        plt.plot(barrel_xx,barrel_yy,color='yellow',linewidth=10)
        plt.plot(bregma_coords[1]-(barrel_xx-bregma_coords[1]),barrel_yy,color='yellow',linewidth=10)      
        plt.plot(motor_xx,motor_yy,color='red',linewidth=10)
        plt.plot(bregma_coords[1]-(motor_xx-bregma_coords[1]),motor_yy,color='red',linewidth=10)
        plt.plot(visual_xx,visual_yy,color='blue',linewidth=10)
        plt.plot(bregma_coords[1]-(visual_xx-bregma_coords[1]),visual_yy,color='blue',linewidth=10)
        plt.plot(rs_xx,rs_yy,color='cyan',linewidth=10)
        plt.plot(bregma_coords[1]-(rs_xx-bregma_coords[1]),rs_yy,color='cyan',linewidth=10)
        plt.plot(acc_xx,acc_yy,color='magenta',linewidth=10)
        plt.plot(bregma_coords[1]-(acc_xx-bregma_coords[1]),acc_yy,color='magenta',linewidth=10)
            
        plt.xlim(0,n_pixels-1)
        plt.ylim(n_pixels-1,0)
        ax.set_title("Labeled areas")
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        plt.show()


def Remove_artifacts(unit, channel, images_processed, file_dir, file_name):
    ''' Tool to manually ablate small regions of imaging area that may be 
        artifacts '''
    
    global coords, images_temp, ax, fig, cid

    print "Manual Mask Mode"
    images_temp = np.array(images_processed).copy()
    
    #Load Generic Mask
    if (os.path.exists(file_dir + 'genericmask.txt')==True):
        print "Loading existing generic mask"
        generic_mask_file = file_dir + file_name + '/' + file_name+ '_genericmask.txt'
        coords_generic = np.loadtxt(generic_mask_file)
        #Ablate generic map
        for i in range(len(images_temp)):
            for j in range(len(coords_generic)):
                images_temp[i][min(255,int(coords_generic[j][0]))][min(255,int(coords_generic[j][1]))]=0
                
    #Load existing specific mask file
    coords=[]
    specific_mask_file = file_dir +'artifactmask.txt'
    if (os.path.exists(specific_mask_file)==True):
        temp_data= np.loadtxt(specific_mask_file)
        for i in range(len(temp_data)):
            coords.append(temp_data[i])
        update_length=len(coords)

        #Ablate specific map
        for i in range(len(images_temp)):
            for j in range(len(coords)):
                for k in range(7):
                    for l in range(7):
                        images_temp[i][min(255,int(coords[j][0])-3+k)][min(255,int(coords[j][1])-3+l)]=0
    else:
        update_length=0

    fig, ax = plt.subplots()
    ax.imshow(images_temp[100])
    cid = fig.canvas.mpl_connect('button_press_event', on_click)
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

    
def Add_tuples(coords, bregma_coords_temp, generic_coords):
    #print coords

    coords = np.array(np.vstack((coords, bregma_coords_temp)), dtype=np.int16)
    coords = np.array(np.vstack((coords, generic_coords)), dtype=np.int16)
    
    #print coords
    #quit()
    
    return coords
    
def Load_areas(unit, channel, n_pixels, file_dir, file_name, images_processed):
        
    #Create set of images for each defined area
    images_areas = []   #Make list of lists to hold images for each area
    coords_save = []    #holds all ablated coords areas
    borders = []
    
    #Load General mask (removes background)
    generic_mask_file = file_dir + 'genericmask.txt'
    generic_coords = np.loadtxt(generic_mask_file)
    
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
                
    #Load Area #1 - Hindlimb
    hindlimb_file = file_dir + 'hindlimb.txt'
    if (os.path.exists(hindlimb_file)==True):
        hindlimb_coords = np.loadtxt(hindlimb_file)     #Load hindlimb coordinates; Mask out everything else
        ipsi_xx = np.append(hindlimb_coords[:,1], hindlimb_coords[0][1])
        ipsi_yy = np.append(hindlimb_coords[:,0], hindlimb_coords[0][0])
        ipsi_coords = np.column_stack((ipsi_yy, ipsi_xx))
        
        contra_coords_xx= bregma_coords[1]-(ipsi_xx-bregma_coords[1])
        contra_coords_yy= ipsi_yy
        contra_coords = np.column_stack((contra_coords_yy, contra_coords_xx))  #Invert these coordinates
        
        ablate_coords1 = Ablate_outside_area(n_pixels, ipsi_coords)
        ablate_coords2 = Ablate_outside_area(n_pixels, contra_coords)
        ablate_coords = np.array([x for x in set(tuple(x) for x in ablate_coords1) & set(tuple(x) for x in ablate_coords2)])

        ablate_coords = Add_tuples(ablate_coords, bregma_coords_temp, generic_coords)
        
        tempmask_indexes=np.zeros((n_pixels,n_pixels))
        for i in range(len(ablate_coords)):
            tempmask_indexes[ablate_coords[i][0]][ablate_coords[i][1]] = True

        a = []
        for i in range(len(images_processed)):
            a.append(np.ma.array(images_processed[i], mask=tempmask_indexes))
        
        images_areas.append(a)
        coords_save.append(ablate_coords)
        borders.append(ipsi_coords)

    #Load Area #2 - Barrel
    barrel_file = file_dir + 'barrel.txt'
    if (os.path.exists(barrel_file)==True):
        barrel_coords = np.loadtxt(barrel_file)     #Load hindlimb coordinates; Mask out everything else
        ipsi_xx = np.append(barrel_coords[:,1], barrel_coords[0][1])
        ipsi_yy = np.append(barrel_coords[:,0], barrel_coords[0][0])
        ipsi_coords = np.column_stack((ipsi_yy, ipsi_xx))
        
        contra_coords_xx= bregma_coords[1]-(ipsi_xx-bregma_coords[1])
        contra_coords_yy= ipsi_yy
        contra_coords = np.column_stack((contra_coords_yy, contra_coords_xx))  #Invert these coordinates

        ablate_coords1 = Ablate_outside_area(n_pixels, ipsi_coords)
        ablate_coords2 = Ablate_outside_area(n_pixels, contra_coords)
        ablate_coords = np.array([x for x in set(tuple(x) for x in ablate_coords1) & set(tuple(x) for x in ablate_coords2)])
        
        ablate_coords = Add_tuples(ablate_coords, bregma_coords_temp, generic_coords)
               
        tempmask_indexes=np.zeros((len(images_processed[0]),len(images_processed[0])))
        for i in range(len(ablate_coords)):
            tempmask_indexes[ablate_coords[i][0]][ablate_coords[i][1]] = True

        a = []
        for i in range(len(images_processed)):
            a.append(np.ma.array(images_processed[i], mask=tempmask_indexes))
        
        images_areas.append(a)
        coords_save.append(ablate_coords)
        borders.append(ipsi_coords)


    #Load Area #3 - Motor
    motor_file = file_dir + 'motor.txt'
    if (os.path.exists(motor_file)==True):
        motor_coords = np.loadtxt(motor_file)     #Load hindlimb coordinates; Mask out everything else
        ipsi_xx = np.append(motor_coords[:,1], motor_coords[0][1])
        ipsi_yy = np.append(motor_coords[:,0], motor_coords[0][0])
        ipsi_coords = np.column_stack((ipsi_yy, ipsi_xx))
        
        contra_coords_xx= bregma_coords[1]-(ipsi_xx-bregma_coords[1])
        contra_coords_yy= ipsi_yy
        contra_coords = np.column_stack((contra_coords_yy, contra_coords_xx))  #Invert these coordinates

        ablate_coords1 = Ablate_outside_area(n_pixels, ipsi_coords)
        ablate_coords2 = Ablate_outside_area(n_pixels, contra_coords)
        ablate_coords = np.array([x for x in set(tuple(x) for x in ablate_coords1) & set(tuple(x) for x in ablate_coords2)])
        
        ablate_coords = Add_tuples(ablate_coords, bregma_coords_temp, generic_coords)
                        
        tempmask_indexes=np.zeros((len(images_processed[0]),len(images_processed[0])))
        for i in range(len(ablate_coords)):
            tempmask_indexes[ablate_coords[i][0]][ablate_coords[i][1]] = True

        a = []
        for i in range(len(images_processed)):
            a.append(np.ma.array(images_processed[i], mask=tempmask_indexes))
        
        images_areas.append(a)
        coords_save.append(ablate_coords)
        borders.append(ipsi_coords)

    #Load Area #4 - Visual
    visual_file = file_dir + 'visual.txt'
    if (os.path.exists(visual_file)==True):
        visual_coords = np.loadtxt(visual_file)     #Load hindlimb coordinates; Mask out everything else
        ipsi_xx = np.append(visual_coords[:,1], visual_coords[0][1])
        ipsi_yy = np.append(visual_coords[:,0], visual_coords[0][0])
        ipsi_coords = np.column_stack((ipsi_yy, ipsi_xx))
        
        contra_coords_xx= bregma_coords[1]-(ipsi_xx-bregma_coords[1])
        contra_coords_yy= ipsi_yy
        contra_coords = np.column_stack((contra_coords_yy, contra_coords_xx))  #Invert these coordinates

        ablate_coords1 = Ablate_outside_area(n_pixels, ipsi_coords)
        ablate_coords2 = Ablate_outside_area(n_pixels, contra_coords)
        ablate_coords = np.array([x for x in set(tuple(x) for x in ablate_coords1) & set(tuple(x) for x in ablate_coords2)])
        
        ablate_coords = Add_tuples(ablate_coords, bregma_coords_temp, generic_coords)
                        
        tempmask_indexes=np.zeros((len(images_processed[0]),len(images_processed[0])))
        for i in range(len(ablate_coords)):
            tempmask_indexes[ablate_coords[i][0]][ablate_coords[i][1]] = True

        a = []
        for i in range(len(images_processed)):
            a.append(np.ma.array(images_processed[i], mask=tempmask_indexes))
        
        images_areas.append(a)
        coords_save.append(ablate_coords)
        borders.append(ipsi_coords)

    #Load Area #5 - Retrosplenial
    rs_file = file_dir + 'rs.txt'
    if (os.path.exists(rs_file)==True):
        rs_coords = np.loadtxt(rs_file)     #Load hindlimb coordinates; Mask out everything else
        ipsi_xx = np.append(rs_coords[:,1], rs_coords[0][1])
        ipsi_yy = np.append(rs_coords[:,0], rs_coords[0][0])
        ipsi_coords = np.column_stack((ipsi_yy, ipsi_xx))
        
        contra_coords_xx= bregma_coords[1]-(ipsi_xx-bregma_coords[1])
        contra_coords_yy= ipsi_yy
        contra_coords = np.column_stack((contra_coords_yy, contra_coords_xx))  #Invert these coordinates

        ablate_coords1 = Ablate_outside_area(n_pixels, ipsi_coords)
        ablate_coords2 = Ablate_outside_area(n_pixels, contra_coords)
        ablate_coords = np.array([x for x in set(tuple(x) for x in ablate_coords1) & set(tuple(x) for x in ablate_coords2)])
        
        ablate_coords = Add_tuples(ablate_coords, bregma_coords_temp, generic_coords)
                        
        tempmask_indexes=np.zeros((len(images_processed[0]),len(images_processed[0])))
        for i in range(len(ablate_coords)):
            tempmask_indexes[ablate_coords[i][0]][ablate_coords[i][1]] = True

        a = []
        for i in range(len(images_processed)):
            a.append(np.ma.array(images_processed[i], mask=tempmask_indexes))
        
        images_areas.append(a)
        coords_save.append(ablate_coords)
        borders.append(ipsi_coords)

    #Load Area #6 - Anterior Cingulate Cortex
    acc_file = file_dir + 'acc.txt'
    if (os.path.exists(acc_file)==True):
        acc_coords = np.loadtxt(acc_file)     #Load hindlimb coordinates; Mask out everything else
        ipsi_xx = np.append(acc_coords[:,1], acc_coords[0][1])
        ipsi_yy = np.append(acc_coords[:,0], acc_coords[0][0])
        ipsi_coords = np.column_stack((ipsi_yy, ipsi_xx))
        
        contra_coords_xx= bregma_coords[1]-(ipsi_xx-bregma_coords[1])
        contra_coords_yy= ipsi_yy
        contra_coords = np.column_stack((contra_coords_yy, contra_coords_xx))  #Invert these coordinates

        ablate_coords1 = Ablate_outside_area(n_pixels, ipsi_coords)
        ablate_coords2 = Ablate_outside_area(n_pixels, contra_coords)
        ablate_coords = np.array([x for x in set(tuple(x) for x in ablate_coords1) & set(tuple(x) for x in ablate_coords2)])
        
        ablate_coords = Add_tuples(ablate_coords, bregma_coords_temp, generic_coords)
                        
        tempmask_indexes=np.zeros((len(images_processed[0]),len(images_processed[0])))
        for i in range(len(ablate_coords)):
            tempmask_indexes[ablate_coords[i][0]][ablate_coords[i][1]] = True

        a = []
        for i in range(len(images_processed)):
            a.append(np.ma.array(images_processed[i], mask=tempmask_indexes))
        
        images_areas.append(a)
        coords_save.append(ablate_coords)
        borders.append(ipsi_coords)

    return images_areas, coords_save, borders, generic_mask_indexes


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

def Search_max_min(images_areas, img_rate, window, coords_save):

    global coords_temp

    area_names = ['hindlimb', 'barrel', 'motor', 'visual', 'retrosplenial', 'acc']

    Max_plot = []
    Min_plot = []
    Max_pixel_value = []
    Min_pixel_value = []
    Images_areas = []

    #Loop over each cortical area
    for a in range(len(images_areas)):
        print "Searching max/min for area # ", area_names[a]
        coords_temp = coords_save[a]

        temp_max_array = []
        temp_min_array = []
        temp_array1 = []
        temp_array2 = []

        for i in range(int(img_rate*2)):     #This searches +/- 1 sec from time = 0 sec
            temp_array1.append(np.array(images_areas[a][int(window*img_rate)-int(img_rate)+i]).copy())
            temp_array2.append(np.array(images_areas[a][int(window*img_rate)-int(img_rate)+i]).copy())

        #Set background values to very large/very low values for comparisons...
        pool = mp.Pool(12)
        temp_max_array.extend(pool.map(Mp_max, temp_array1))
        temp_min_array.extend(pool.map(Mp_min, temp_array2))
        pool.close()

        temp_max_array = np.array(temp_max_array)
        max_index = np.unravel_index(np.argmax(temp_max_array), temp_max_array.shape)
        max_pixel_value = images_areas[a][int(window*img_rate)-int(img_rate)+max_index[0]][max_index[1]][max_index[2]]

        temp_min_array = np.array(temp_min_array)
        min_index = np.unravel_index(np.argmin(temp_min_array), temp_min_array.shape)
        min_pixel_value = images_areas[a][int(window*img_rate)-int(img_rate)+min_index[0]][min_index[1]][min_index[2]]
        
        max_index = max_index[1:] #Reduce back to 2D arrays
        min_index = min_index[1:]
       
        max_plot = []
        min_plot = []
        
        for i in range(len(images_areas[a])):
            #Time course curve data - before zeroing out the pixels
            max_plot.append(images_areas[a][i][max_index[0]][max_index[1]])
            min_plot.append(images_areas[a][i][min_index[0]][min_index[1]])

            #Zero out pixels - max index
            for p in range(3):
                for r in range(3):
                    if (max_index[0]-1+p)<256 and (max_index[0]-1+p)>0 and (max_index[1]-1+r)<256 and (max_index[1]-1+r)>0:
                        images_areas[a][i][max_index[0]-1+p][max_index[1]-1+r]=max_pixel_value

            #Zero out pixels - min index
            for p in range(3):
                for r in range(3):
                    images_areas[a][i][min_index[0]-1+p][min_index[1]-1+r]=max_pixel_value
          
        Max_plot.append(max_plot)
        Min_plot.append(min_plot)
        Max_pixel_value.append(max_pixel_value)
        Min_pixel_value.append(min_pixel_value)
        Images_areas.append(images_areas[a])
    
    return Images_areas, Max_plot, Min_plot, Max_pixel_value, Min_pixel_value

def Save_time_course(unit, channel, spikes, max_plot, min_plot, window, len_frame, file_dir, file_name):
    
    colors=['blue','red','green']

    max_plot=np.array(max_plot)
    min_plot=np.array(min_plot)

    xx = np.arange(-window, window, len_frame)
    xx = xx[0: len(max_plot)]

    fig, ax=plt.subplots(nrows=1,ncols=1)
    ax.plot(xx, max_plot, color='black', linewidth=2)
    ax.plot(xx, min_plot, color='blue', linewidth=2)
    ax.plot([-3,3],[0,0], color='black')
    ax.plot([0,0],[min(min_plot),max(max_plot)], color='black')

    plt.suptitle(file_name + "    Max/Min pixel - time course \n Unit: " +str(unit).zfill(2) + " Channel: " + str(channel).zfill(2)+
    " \n No. of trigger spikes: "+ str(len(spikes)))
    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    fig.savefig(file_dir + file_name + '/' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'.png', fontsize = 20)
    plt.close(fig)
    
    temp_array = []
    temp_array.append(min_plot)
    temp_array.append(max_plot)
    temp_array = np.array(temp_array, dtype=np.float32)
    np.save(file_dir + file_name + '/time_course_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2), temp_array)

# Draw frames and write to the pipe
def Par_anim((args)):
    
    global plt_string, unt, spks, _r, generic_coords, f_dir, f_name, Min_pixel_v, Max_pixel_v, Max_p, Min_p, win, colorz, bregma_c, borderz, n_pix, len_fr

    index, Images_areas = args
    print "Index: ", index

    fig = plt.figure()
    gs = gridspec.GridSpec(50,110)
    axes = []
    axes.append(fig.add_subplot(gs[:,:50]))
    axes.append(fig.add_subplot(gs[10:18,60:80]))
    axes.append(fig.add_subplot(gs[21:29,60:80]))
    axes.append(fig.add_subplot(gs[32:41,60:80]))
    axes.append(fig.add_subplot(gs[10:18,90:110]))
    axes.append(fig.add_subplot(gs[21:29,90:110]))
    axes.append(fig.add_subplot(gs[32:41,90:110]))

    xx = np.arange(-win, win, len_fr)
    xx = xx[0: len(Max_p[0])]
 
    masked_data = np.ma.array(Images_areas[0], mask=generic_coords)
    axes[0].imshow(masked_data, vmin=min(Min_pixel_v), vmax=max(Max_pixel_v), origin='lower') #cmap=my_cmap, clim=[0.9, 1]) #, cmap=cmap, interpolation='nearest')

    for a in range(len(Images_areas)):
        axes[0].plot(borderz[a][:,1],borderz[a][:,0],color=colorz[a],linewidth=3)
        axes[0].plot(bregma_c[1]-(borderz[a][:,1]-bregma_c[1]),borderz[a][:,0],color=colorz[a],linewidth=3)
        axes[0].set_xlim(0,n_pix-1)
        axes[0].set_ylim(n_pix-1,0)
    
        axes[a+1].plot(xx[0:index], Max_p[a][0:index], color=colorz[a], linewidth=2)
        axes[a+1].plot(xx[0:index], Min_p[a][0:index], ':', color=colorz[a], linewidth=2)
        axes[a+1].plot([0,0], [-1,+1], color='black', linewidth = .5, alpha=0.5)
        axes[a+1].plot([-win,+win], [0,0], color='black', linewidth = .5, alpha=0.5)
        axes[a+1].set_ylim(-0.03, 0.03)
        axes[a+1].tick_params(axis='both', which='major', labelsize=10)
        
    plt.suptitle(f_name + " Unit: "+str(unt) + " No. of trigger spikes: "+ str(len(spks))+ '\n'+ "Spikes plotted: " + plt_string)
    axes[0].set_title("Time: %.2f secs." % (float(index)/img_r-win))
    #ax1.text(3, 8, "Time : "+str(float(time)/img_r-win), style='italic') #bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})

    fig.canvas.draw()
    fig.savefig(f_dir + f_name+'/figs_'+str(index).zfill(3)+'.png')
    plt.close(fig)
    #string = fig.canvas.tostring_argb()
    #return string
    
    
def Animate_images(unit, channel, window, len_frame, Max_plot, Min_plot, Images_areas, file_dir, file_name, 
    Max_pixel_value, Min_pixel_value, n_pixels, borders, generic_mask_indexes, img_rate, spikes, plot_string):

    global plt_string, unt, spks, img_r, win, generic_coords, f_dir, f_name, index_array, Min_pixel_v, Max_pixel_v, Max_p, Min_p, win, colorz, bregma_c, borderz, n_pix, len_fr

    print "Generating Animation"
    
    #Load Bregma coords for computing contralateral areas
    bregma_mask_file = file_dir + 'bregmamask.txt'
    bregma_coords = np.loadtxt(bregma_mask_file)
    
    generic_coords = generic_mask_indexes

    colors = ['green', 'brown', 'red', 'blue', 'cyan', 'magenta']

    Min_pixel_v = Min_pixel_value
    Max_pixel_v = Max_pixel_value
    Max_p = Max_plot
    Min_p = Min_plot
    win = window
    colorz = colors
    bregma_c = bregma_coords
    borderz = borders
    n_pix = n_pixels
    len_fr = len_frame
    index_array=0
    f_dir = file_dir
    f_name = file_name
    win = window
    img_r = img_rate
    spks = spikes
    unt = unit
    plt_string = plot_string

    strings=[]
    Images_areas = np.array(Images_areas)
    Images_areas = np.swapaxes(Images_areas, 0,1)
    #Images_areas = Images_areas[0:5]

    indices = np.arange(len(Images_areas))
    
    pool = mp.Pool(12)
    pool.map(Par_anim, zip(indices, Images_areas))
    pool.close()
    
    print "Excuting command", 'ffmpeg -f image2 -r 15 -i ' + file_dir+file_name+'/figs_%03d.png -vcodec libx264 -y '+file_dir+file_name+'/'+file_name+'_unit_'+str(unit)+'_'+plot_string+'.mp4'

    os.system('ffmpeg -f image2 -r 15 -i ' + file_dir+file_name+'/figs_%03d.png -vcodec libx264 -y '+file_dir+file_name+'/'+file_name+
    '_unit_'+str(unit)+'_'+plot_string+'.mp4')
    os.system('rm '+file_dir+file_name+'/figs*')

    print "Finished unit: ", unit

    #if False:

        #Writer = animation.writers['ffmpeg']
        #writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

        ##setup figure
        #fig = plt.figure()
        #gs = gridspec.GridSpec(50,110)
        #axes = []
        #axes.append(fig.add_subplot(gs[:,:50]))
        #axes.append(fig.add_subplot(gs[10:18,60:80]))
        #axes.append(fig.add_subplot(gs[21:29,60:80]))
        #axes.append(fig.add_subplot(gs[32:41,60:80]))
        #axes.append(fig.add_subplot(gs[10:18,90:110]))
        #axes.append(fig.add_subplot(gs[21:29,90:110]))
        #axes.append(fig.add_subplot(gs[32:41,90:110]))

        ##set up list of images for animation
        #ims=[]
        #xx = np.arange(-window, window, len_frame)
        #xx = xx[0: len(Max_plot[0])]
    
        #for time in range(1): #len(Images_areas[0])):
            #print "Time: ", time
            #ims_temp = []
            ##Add area specific plots
            #for a in range(len(Images_areas)):
                #im_00 = axes[0].imshow(Images_areas[a][time], vmin=min(Min_pixel_value), vmax=max(Max_pixel_value), interpolation='none')
                #axes[0].plot(borders[a][:,1],borders[a][:,0],color=colors[a],linewidth=3)
                #axes[0].plot(bregma_coords[1]-(borders[a][:,1]-bregma_coords[1]),borders[a][:,0],color=colors[a],linewidth=3)
                #axes[0].set_xlim(0,n_pixels-1)
                #axes[0].set_ylim(n_pixels-1,0)
            
                #im_01, = axes[a+1].plot(xx[0:time], Max_plot[a][0:time], color=colors[a], linewidth=2)
                #im_02, = axes[a+1].plot(xx[0:time], Min_plot[a][0:time], ':', color=colors[a], linewidth=2)
                #im_03, = axes[a+1].plot([0,0], [-1,+1], color='black', linewidth = .5, alpha=0.5)
                #im_04, = axes[a+1].plot([-window,+window], [0,0], color='black', linewidth = .5, alpha=0.5)
                
                #axes[a+1].set_ylim(-0.03, 0.03)

                #ims_temp.extend([im_00, im_01, im_02, im_03, im_04])
            
            #ims.append(ims_temp)


        #ani = anim.ArtistAnimation(fig, ims, interval=50, blit=True)
        #print "Saving Animation"

        #ani.save(file_dir + file_name + '/' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'.mp4', writer=writer)


    ##Subprocess option opens PIPE and feeds 'strings' to the ffmpeg compiler via command line
    #if False: 
        ## Open an ffmpeg process
        #canvas_width = 640
        #canvas_height = 480
        #outf = file_dir + file_name + '/' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'.mp4'
        #cmdstring = ('ffmpeg', 
                     #'-y', '-r',  '15', # overwrite, frame rate
                     #'-s', '%dx%d' % (canvas_width, canvas_height), # size of image string
                     #'-pix_fmt', 'argb', # format
                     #'-f', 'rawvideo',  '-i', '-', # tell ffmpeg to expect raw video from the pipe
                     #'-vcodec', 'mpeg4',
                      #outf) # output file
              
                     
        #import subprocess
        #from subprocess import PIPE
        #p = subprocess.Popen(cmdstring, stdin=PIPE, stdout=PIPE)
        #pool = mp.Pool(12)
        #strings.extend(pool.map(Par_anim, zip(indices, Images_areas)))
        
        #print shape(strings)
        #for i in range(len(strings)):
            #p.stdin.write(strings[i])
    

   # pool = mp.Pool(12)
   # temp_min_array.extend(pool.map(Mp_min, temp_array2))
   # pool.close()

    #def Mp_min((temp_array2)):
        #global coords_temp

        #for j in range(len(coords_temp)):
            #temp_array2[coords_temp[j][0]][coords_temp[j][1]] = 1000

        #return temp_array2
    
    # Finish up
    #p.communicate()



#def konnerth_baseline(frames): 
    
    #temp_array=np.zeros((128, 128), dtype=np.float32)+255
    #temp_sum = 1E10
    #if True:
        #for i in range(len(frames)/2-25):
            #temp_baseline = np.sum(np.mean(frames[i:i+25], axis=0))
            ##print temp_baseline
            #if temp_baseline < temp_sum: 
                #temp_sum = temp_baseline
                #temp_array = np.mean(frames[i:i+25], axis=0)

    ##if True: temp_array = np.mean(frames[0:len(frames)/2],axis=0)

    #return temp_array

#def konnerth_filtering(frames): 
    #temp_frames = np.zeros((int(window*img_rate)*2, n_pixels, n_pixels), dtype=np.float32)
    
    #for i in range(int(window*img_rate)*2-7):
        #temp_exp = 0
        #temp_sum = np.zeros((n_pixels, n_pixels), dtype=np.float32)
        #for t in range(5):
            #temp_sum += np.exp(-5*t) * frames[i+t]
            #temp_exp += np.exp(-5*t)
        #temp_frames[i] = temp_sum/temp_exp

    #return temp_frames
