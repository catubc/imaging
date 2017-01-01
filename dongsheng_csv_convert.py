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

locations = ['hindlimb_left','hindlimb_right','forelimb_left','forelimb_right','barrel_left','barrel_right','motor_left',
'motor_right','visual_left','visual_right','retrosplenial_left','retrosplenial_right','acc_left','acc_right','allcortex']

main_dir = '/media/cat/12TB/in_vivo/tim/'

file_dirs = []
file_names = []

#****************************************************************************************************************************
if False:
    file_dirs.append('/media/cat/12TB/in_vivo/tim/2015-11-18/') #*************************************
    file_names.append([
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

    file_dirs.append('/media/cat/12TB/in_vivo/tim/2015-11-27/') #*************************************
    file_names.append([
    '2015-11-27-4-10electrodein-iso0', 
    '2015-11-27-5-10electrodein-iso0',
    '2015-11-27-10-16electrodein-iso0',
    '2015-11-27-11-16electrodein-iso0',
    '2015-11-27-16-deep-iso0',
    ])

    file_dirs.append('/media/cat/12TB/in_vivo/tim/2015-12-1/') #*************************************
    file_names.append([
    '2015-12-1-5-10electrodeiniso0',
    '2015-12-1-6-10electrodeiniso0',
    '2015-12-1-11-allelectrodeiniso0', 
    '2015-12-1-12-allelectrodeiniso0',
    '2015-12-1-18-5electrodeinthalamus-iso0',   #SINGLE
    '2015-12-1-22-allelectrodeinthalamus-iso0',
    '2015-12-1-23-allelectrodeinthalamus-iso0', 
    ##'2015-12-1-24-allelectrodeinthalamus-iso0'
    ])


    file_dirs.append('/media/cat/12TB/in_vivo/tim/2015-12-2/') #*************************************
    file_names.append([
    '2015-12-2-6-10electrodeincortex-iso0',
    '2015-12-2-7-10electrodeincortex-iso0',
    '2015-12-2-14-allelectrodeinthalamus-is0',
    '2015-12-2-15-allelectrodeinthalamus-is0',
    ])


    file_dirs.append('/media/cat/12TB/in_vivo/tim/2015-12-11/') #*************************************
    file_names.append([
    '2015-12-11-5-10electincortex-iso0',
    '2015-12-11-7-10electincortex-iso1',
    '2015-12-11-15-allelectinthalamus-iso0',
    '2015-12-11-16-allelectinthalamus-iso0',
    '2015-12-11-17-allelectinthalamus-iso0',
    ])

    file_dirs.append('/media/cat/12TB/in_vivo/tim/2015-12-16/') #*************************************
    file_names.append([
    '2015-12-16-4-11electrodeincortex-iso0',
    '2015-12-16-5-11electrodeincortex-iso0',
    '2015-12-16-13-allelectrodeinZI-iso0',
    '2015-12-16-14-allelectrodeinZI-iso0',
    ])

    file_dirs.append('/media/cat/12TB/in_vivo/tim/2016-1-11/') #*************************************
    file_names.append([
    '2016-1-11-5-10electrodeincortex-iso0',
    '2016-1-11-6-10electrodeincortex-iso0',
    '2016-1-11-14-13electrodeinthalamus-iso0',
    '2016-1-11-15-13electrodeinthalamus-iso0',
    ])

    file_dirs.append('/media/cat/12TB/in_vivo/tim/2016-1-14/') #*************************************
    file_names.append([
    '2016-1-14-5-10electrodein cortex-iso0',
    '2016-1-14-6-10electrodein cortex-iso0',
    '2016-1-14-7-10electrodein cortex-iso0',
    '2016-1-14-17-allelectrodeinthalamus-iso0'
    ])


#file_dirs.append('/media/cat/12TB/in_vivo/tim/2015-5-6/')  #*************************************
#file_names.append([
#'2015-5-6-1',       #cortex
#'2015-5-6-2',       #cortex
#'2015-5-6-3',      #subcortex
#'2015-5-6-4',      #subcortex
#'2015-5-6-5',      #subcortex
#'2015-5-6-6'       #subcortex
#])

##file_dirs.append('/media/cat/12TB/in_vivo/tim/2015-7-22/')  #*************************************
##file_names.append([
##'2015-7-22-1',     #cortex
##'2015-7-22-2',     #cortex
##'2015-7-22-3',     #cortex
##'2015-7-22-4',     #subcortex
##'2015-7-22-5',     #subcortex
##'2015-7-22-6',     #subcortex
##'2015-7-22-7',     #subcortex
##'2015-7-22-8'      #subcortex
##])

#file_dirs.append('/media/cat/12TB/in_vivo/tim/2015-7-23/')  #*************************************
#file_names.append([
#'2015-7-23-1',     #cortex
#'2015-7-23-2',     #subcortex
#'2015-7-23-3',     #subcortex  
#'2015-7-23-16'     #subcortex
#])


##file_dirs.append('/media/cat/12TB/in_vivo/tim/2015-7-30/')  #*************************************
##file_names.append([
##'2015-7-30-4',
##'2015-7-30-5', 
##'2015-7-30-6', 
##'2015-7-30-7'
##])

file_dirs.append('/media/cat/12TB/in_vivo/tim/2015-11-18/') #*************************************
file_names.append([
'2015-11-18-1-9electrodein', 
'2015-11-18-2-9electrodein', 
'2015-11-18-3-9electrodein',
'2015-11-18-4-9electrodein-iso0.5',
'2015-11-18-5-9electrodein-iso0.5',
'2015-11-18-6-9electrodein-iso0.5',
])


file_dirs.append('/media/cat/12TB/in_vivo/tim/2015-11-27/') #*************************************
file_names.append([
'2015-11-27-1-10electrodein-iso1.0',
'2015-11-27-2-10electrodein-iso1', 
'2015-11-27-7-16electrodein-iso1',
'2015-11-27-8-16electrodein-iso1',
'2015-11-27-13-deep-iso1',
'2015-11-27-14-deep-iso1',
])

file_dirs.append('/media/cat/12TB/in_vivo/tim/2015-12-1/') #*************************************
file_names.append([
'2015-12-1-1-10electrodeiniso1',
'2015-12-1-2-10electrodeiniso1',
'2015-12-1-3-10electrodeiniso1',
'2015-12-1-8-allelectrodeiniso0.8', 
'2015-12-1-9-allelectrodeiniso0.8', 
'2015-12-1-14-5electrodeinthalamus-iso0.8',
'2015-12-1-15-5electrodeinthalamus-iso0.8',
'2015-12-1-19-allelectrodeinthalamus-iso0.8',  #SINGLE
'2015-12-1-20-allelectrodeinthalamus-iso0.8',
])

file_dirs.append('/media/cat/12TB/in_vivo/tim/2015-12-2/') #*************************************
file_names.append([
'2015-12-2-2-10electrodeincortex-iso1', 
'2015-12-2-3-10electrodeincortex-iso1',
'2015-12-2-4-10electrodeincortex-iso1',
'2015-12-2-9-allelectrodeincortex-iso1',  #SINGLE
'2015-12-2-11-allelectrodeinthalamus-iso1',
'2015-12-2-12-allelectrodeinthalamus-iso1',
])


file_dirs.append('/media/cat/12TB/in_vivo/tim/2015-12-11/') #*************************************
file_names.append([
'2015-12-11-1-10electincortex-iso1.0',
'2015-12-11-2-10electincortex-iso1',
'2015-12-11-3-10electincortex-iso1',
'2015-12-11-7-10electincortex-iso1',
'2015-12-11-8-10electincortex-iso1',
'2015-12-11-9-allelectincortex',
'2015-12-11-10-5electinthalamusorcpu',
'2015-12-11-11-allelectinthalamus-iso1',
'2015-12-11-12-allelectinthalamus-iso1',
'2015-12-11-13-allelectinthalamus-iso1',
'2015-12-11-18-allelectinthalamus-iso1'
])

file_dirs.append('/media/cat/12TB/in_vivo/tim/2015-12-16/') #*************************************
file_names.append([
'2015-12-16-1-11electrodeincortex-iso1',
'2015-12-16-2-11electrodeincortex-iso1',
'2015-12-16-7-11electrodeincortex-iso1',
'2015-12-16-8-4electrodeincpu-iso1',
'2015-12-16-9-allelectrodeinZI-iso1',
'2015-12-16-10-allelectrodeinZI-iso1',
'2015-12-16-11-allelectrodeinZI-iso1',
'2015-12-16-16-allelectrodeinZI-iso1',
'2015-12-16-17-allelectrodeinZI-iso1'
])


file_dirs.append('/media/cat/12TB/in_vivo/tim/2016-1-11/') #*************************************
file_names.append([
'2016-1-11-1-10electrodeincortex-iso1.2',
'2016-1-11-2-10electrodeincortex-iso1',
'2016-1-11-3-10electrodeincortex-iso1',
'2016-1-11-8-10electrodeincortex-iso1',
'2016-1-11-9-allelectrodeincortex-iso1',
'2016-1-11-10-13electrodeinthalamus-iso1',
'2016-1-11-11-13electrodeinthalamus-iso1',
'2016-1-11-12-13electrodeinthalamus-iso1',
'2016-1-11-17-13electrodeinthalamus-iso1',
'2016-1-11-18-13electrodeinthalamus-iso1'
])



file_dirs.append('/media/cat/12TB/in_vivo/tim/2016-1-14/') #*************************************
file_names.append([
'2016-1-14-1-10electrodein cortex-iso1',
'2016-1-14-2-10electrodein cortex-iso1',
'2016-1-14-3-10electrodein cortex-iso1',
'2016-1-14-9-10electrodein cortex-iso1',
'2016-1-14-10-allelectrodein cortex-iso1',
'2016-1-14-11-allelectrodeinthalamus-iso1',
'2016-1-14-12-allelectrodeinthalamus-iso1',
'2016-1-14-13-allelectrodeinthalamus-iso1',
])


#Quality control parameters

minDF_response = 0.005
min_spikes = 200

total_cell_count = 0

dir_counter = 0
for file_dir in file_dirs:
    print file_dir
    for file_name in file_names[dir_counter]:
        print file_name
        
        #Load quality control file; exclude units that did not make QC!
        if False:
            temp_fname = file_dir.replace(main_dir, '')
            dir_file_name = main_dir+"dongsheng_quality_control/"+temp_fname+file_name+'/'
            files = os.listdir(dir_file_name)
            qc = np.zeros(100,dtype=np.int32)
            for f1 in files:
                #print f1, f1[-6:-4]
                qc[int(f1[-6:-4])] = 1  #S


        #Load units using .csv file name
        files = os.listdir(file_dir+file_name)
        temp_names = []
        for file_ in files:
            if ("unit_" in file_):
                if ('map' in file_):
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
            if False: #SKIP QUALITY CONTROL - as per Dongsheng request
                if qc[temp_unit_id]==0: continue

            #Load good units
            units.append(temp_unit_id)
            channels.append(int(file_[16:18]))
            ptps.append(int(file_[23:26]))

        print "No. units: ", len(units)
        
        #Load .tsf raw file to show raw ECP
        tsf_dir = file_dir + file_name+ '/'
        
        tsf_name = tsf_dir + file_name + '_hp.tsf'
        f = open(tsf_name, "rb")
        print "Loading "+ tsf_name 

        tsf = Tsf_file(f, tsf_dir)  #Auto load tsf file attributes: n_electrodes, ec_traces, SampleFreqeuncy and others
        tsf.sim_dir = tsf_dir
        tsf.tsf_name = tsf_name
        
        
        
        #Load beginning and end of recording data
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

            #Load rec_index
            rec_file = file_dir + file_name+ '/rec_index.txt'
            rec_index = np.loadtxt(rec_file)
            #print rec_index
            #print start_array, end_array
            start_array = [start_array[int(rec_index)-1]]; end_array= [end_array[int(rec_index)-1]]

        else:
            epoch_file = file_dir + file_name+ '/epochs.txt'
            data = np.loadtxt(epoch_file)
            start_array = data[:,0]
            end_array = data[:,1]
            
            #Load rec_index
            rec_file = file_dir + file_name+ '/rec_index.txt'
            rec_index = np.loadtxt(rec_file)
            start_array = [start_array[int(rec_index)-1]]; end_array= [end_array[int(rec_index)-1]]

        
        #Load sorted spikes for Dongsheng 
        if True:
            for i in range(0, len(units),1):

                print "*Processing ", file_name, " unit: ", i+1, " of ", len(units), " unit name: ", units[i]
                unit=units[i]
                channel=channels[i]
                ptp=ptps[i]

                if True:
                    in_file = file_dir+file_name+'/'+'unit_'+str(unit).zfill(2)+ '_channel_' + str(channel).zfill(2) + '_ptp_'+str(ptp).zfill(3)
                    spikes = np.loadtxt(in_file+'.csv')
                    
                    spikes = spikes[np.where((spikes>start_array) & (spikes<end_array))[0]]

                    #Quality control - # spikes
                    if len(spikes)<min_spikes: 
                        print "Too few spikes..."
                        continue

                    #Quality control - min DF/F
                    plot_string='all'
                    window = 3
                    tc_name = glob.glob(file_dir + file_name+'/time_course_data_'+file_name+'_'+plot_string+'_'+str(window)+'sec_window_unit'+str(unit).zfill(2)+'*')

                    
                    file_ = tc_name[0]
                    if ("time_course_data" in file_) and ("15sec" not in file_) and ("all_" in file_):

                        stmtd=[]
                        with open(tc_name[0], "r") as f:
                            data = csv.reader(f)
                            data.next()
                            data.next()
                            counter=2
                            for row in data: #skip rows until get to pixel locations
                                stmtd.append([float(i) for i in row])
                                counter+=1
                                if counter>31: break

                        #temp_DF_over_F = max(max(stmtd[29]),abs(min(stmtd[28])))   #Use all-cortex max
                        temp_DF_over_F = max(max(stmtd[9]),abs(min(stmtd[8])))      #Use only barrel left max/mins
                        if temp_DF_over_F<minDF_response: 
                            print "DF/F too low..."
                            continue

                        spikes = np.insert(spikes,0,start_array)
                        spikes = np.insert(spikes,len(spikes),end_array)
                        
                        if True:  #Save new spike files for dongsheng
                            home_dir = '/home/cat/Dropbox/tim/rasters/'
                            exp_name = file_dir.replace('/media/cat/12TB/in_vivo/tim/', '')
                            out_file = home_dir + exp_name+file_name+'/'
                            temp2 = in_file.replace('/media/cat/12TB/in_vivo/tim/'+exp_name+file_name+'/','')
                            out_file = out_file + temp2 +"_ds_.csv"
                            np.savetxt(out_file, spikes)

                        total_cell_count+=1
            print "Done saving dongsheng files..."

    print "Total_cell_count update: ", total_cell_count
    dir_counter+=1
    



if False:

        #Load original .tif file and determine # of frames
        if False:
            img = Image.open(file_dir + file_name + '/' + file_name+'.tif')
            print "Total img frames: ", img.n_frames

            #Load processed shifted images
            shifted_image = np.load(file_dir + file_name+'/'+file_name+'_images_shifted.npy')

            #Remove baseline - global activity regression
            if True:
                print "Removing baseline over all data"
                images_temp = shifted_image.copy()
                #baseline = np.mean(shifted_image, axis=0)
                baseline = np.mean(images_temp, axis=0, dtype=np.float32)
                shifted_image = (shifted_image - baseline)/baseline    
                
        #Load starting image frames to match with ephys time line
        temp_name = glob.glob(file_dir + file_name+'/'+file_name+'_images_start_' + '*')
        end_loc = temp_name[0].index('end')
        last_frame = int(temp_name[0][end_loc+4:])
        
        temp_name[0]=temp_name[0][:end_loc-1]
        start_loc = temp_name[0].index('start')
        first_frame = int(temp_name[0][start_loc+6:])
        
        print "Start/end of imaging: ", start_array[0], end_array[0]
        print "Duration of all imaging: ", end_array[0] - start_array[0], " secs."
        print "Total frames (original .tif): ", img.n_frames
        img_rate = img.n_frames/(end_array[0] - start_array[0])
        print "Imaging rate: ", img_rate
        
        print "Begin/end of brain imaging : ", first_frame, last_frame
        img_begin = first_frame/img_rate
        img_end = last_frame/img_rate
        print "Brain imaging begin and end: ", img_begin, img_end
        
        print "Ephys beging and end: ", start_array, end_array
                
        #**************** LOOP OVER ALL UNITS IN EXPERIMENT **************************
        
        for i in range(0, len(units),1):
            gs = gridspec.GridSpec(10,20)

            print "*Processing ", file_name, " unit: ", i+1, " of ", len(units), " unit name: ", units[i]
            unit=units[i]
            channel=channels[i]
            ptp=ptps[i]

            print "Unit: ", unit
            #if unit in cell_list:
                #pass
            #else: 
                #continue
            
            #Load sorte dspikes for Dongsheng 
            if True:
                in_file = file_dir+file_name+'/'+'unit_'+str(unit).zfill(2)+ '_channel_' + str(channel).zfill(2) + '_ptp_'+str(ptp).zfill(3)
                spikes = np.loadtxt(in_file+'.csv')
                
                spikes = spikes[np.where((spikes>start_array) & (spikes<end_array))[0]]

                spikes = np.insert(spikes,0,start_array)
                spikes = np.insert(spikes,len(spikes),end_array)
                


            #Load [Ca] and plot vs. rasters vs. raw extracellular data 
            if True:
                #Load all max/min pixel locations
                plot_string='all'
                window = 3
                tc_name = glob.glob(file_dir + file_name+'/time_course_data_'+file_name+'_'+plot_string+'_'+str(window)+'sec_window_unit'+str(unit).zfill(2)+'*')
                
                pixel_locations=[]
                stmtd=[]
                with open(tc_name[0], "r") as f:
                    data = csv.reader(f)
                    data.next()
                    data.next()
                    counter=2
                    for row in data: #skip rows until get to pixel locations
                        stmtd.append([float(i) for i in row])
                        counter+=1
                        if counter>31: break
                    for row in data:
                        pixel_locations.append([int(i) for i in row])

                tc_x = np.linspace(-3,3,len(stmtd[28]))
                ax = plt.subplot(gs[1:5,15:20])
                plt.plot(tc_x, np.array(stmtd[28])*1E2, linewidth = 3, color='blue')
                plt.plot([0,0],[-.05,.10])
                plt.plot([-3,3], [0,0])
                plt.ylim(-5,10)
                
                ax = plt.subplot(gs[6:10,15:20])
                plt.plot(tc_x,np.array(stmtd[29])*1E2, linewidth = 3, color='blue')
                plt.ylim(-5,10)
                plt.plot([-3,3], [0,0])
                plt.plot([0,0],[-.05,.10])
                
                time_courses=[]
                for k in range(30):
                    time_courses.append([])

                xx = np.linspace(img_begin+start_array[0],img_end+start_array[0],len(shifted_image))
                

                for location in [28,29]:
                    print "Location: ", location
                    
                    if (location %2)==0:
                        ax = plt.subplot(gs[0:5,0:14])
                    else:
                        ax = plt.subplot(gs[6:10,0:14])

                    for frame in range(len(shifted_image)):
                        pixel_val = shifted_image[frame][pixel_locations[location][0]][pixel_locations[location][1]]
                        #print pixel_val
                        time_courses[location].append(pixel_val)
                        
                    #Plot [Ca] traces
                    plt.plot(xx, np.array(time_courses[location]), linewidth = 3, color='red')

                    ax.get_yaxis().set_ticks([])
                    if (location==0): plt.title("Min Values")
                    if (location==1): plt.title("Max Values")
                    if (location % 2)==0: 
                        h = plt.ylabel(locations[location/2], labelpad=60)
                        h.set_rotation(0)
    
                    #Plot spike rasters
                    for s in spikes:
                        plt.vlines(s, -.175,-.125, color='blue')
                
                    #Plot raw highpass signal
                    start_hp = int(start_array[0]*tsf.SampleFrequency)
                    end_hp = int(end_array[0]*tsf.SampleFrequency)
                    
                    xxx = np.linspace(start_array[0],end_array[0], (end_array[0]-start_array[0])*tsf.SampleFrequency)
                    yyy = tsf.ec_traces[channel-1][start_hp:end_hp]
                    xxx = xxx[:min(len(xxx),len(yyy))]
                    yyy = yyy[:min(len(xxx),len(yyy))]
                    
                    plt.plot(xxx, np.float32(yyy)/tsf.SampleFrequency/2-.3, color='black', alpha=.4)
                    
                    plt.xlabel("Seconds")

                #plt.tight_layout()
                plt.suptitle(file_name+ "   Unit: "+str(unit)+"   #spikes: " + str(len(spikes)))
                plt.show()
            #quit()

