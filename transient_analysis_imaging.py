#from tsf_ptcs_classes import *
#from sequential_firing import *
#from opengl_pyqt_classes import *
#from sta_utils import *
import filter
import os

import numpy as np
import matplotlib.pyplot as plt
import struct
import string
import matplotlib.animation as animation
from scipy import ndimage
from PIL import Image

np.set_printoptions(formatter={'float': '{: 0.6f}'.format})

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
    

def load_leverpositions(lever_position_filename, pull_filename, images_raw):

    import re

    #Load code file
    text_file = open(pull_filename+'.txt', "r")
    lines = text_file.read().splitlines()
    code_text = []
    for line in lines:
        code_text.append(re.split(r'\t+',line))
    del code_text[0:2]
    
    #Load lever position file time stamps
    text_file = open(lever_position_filename+'.txt', "r")
    lines = text_file.read().splitlines()
    lever_text = []
    for line in lines:
        lever_text.append(re.split(r'\t+',line))

    lever_text = np.array(lever_text)
    counter=0
    for k in range(len(lever_text)):
        if lever_text[k][0][:3]=='201':  
            print lever_text[k], code_text[counter]
            counter+=1
    quit()

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
    
#*********************** DATA FILES *******************
save_animation = False
spontaneous = False
stimulus_evoked = False
yuki = True

#Spontaneous 
if spontaneous:
    n_frames = 5000
    offset_frames = 1000

    work_dir = '/media/cat/12TB/in_vivo/tim/allen/20150303/Spontaneous/'
    file_name = '02_sponVSD_150Hz_10000fr_iso0.50%'
    file_names = [file_name]
    file_names_control = []

#Stimulus evoked 
if stimulus_evoked:
    n_frames = 108
    offset_frames = 0

    work_dir = '/media/cat/12TB/in_vivo/tim/allen/20150303/Evoked/L-HL_stim_1ms_1mA/'
    #work_dir = '/media/cat/12TB/in_vivo/tim/allen/20150303/Evoked/L-FL stim_1ms 1mA/'
    #work_dir = '/media/cat/12TB/in_vivo/tim/allen/20150303/Evoked/L-FL stim_1ms 2mA/'
    #work_dir = '/media/cat/12TB/in_vivo/tim/allen/20150303/Evoked/L-whisker stim_1ms 5V/'
    #work_dir = '/media/cat/12TB/in_vivo/tim/allen/20150303/Evoked/L-visual stim_1ms/'
    #work_dir = '/media/cat/12TB/in_vivo/tim/allen/20150303/Evoked/L-aud stim_1ms/'

    file_names = ['hi1', 'hi2', 'hi3', 'hi4', 'hi5','lo1', 'lo2', 'lo3', 'lo4', 'lo5']
    file_names_control = ['no1', 'no2', 'no3', 'no4', 'no5','no1', 'no2', 'no3', 'no4', 'no5']


#Yuki's leaver pull data
if yuki: 
    n_frames = 100000
    offset_frames = 0

    #********** I1 Mouse ************
    #work_dir = '/media/cat/12TB/in_vivo/tim/yuki/I1/'
    
    #file_names = ['I1pm_Apr15_30Hz']
    #pull_filename = work_dir + 'I1_4-15-2015_PM'
    
    #file_names = ['I1am_Apr16_30Hz']
    #pull_filename = work_dir + 'I1_4-16-2015_AM'
    
    #file_names = ['I1am_Apr22_30Hz']
    #pull_filename = work_dir + 'I1_4-22-2015_AM'
    
    #file_names = ['I1am_Apr29_30Hz']
    #pull_filename = work_dir +  'I1_4-29-2015_AM'
    
    #file_names = ['I1am_Apr28_30Hz']
    #pull_filename = work_dir +  'I1_4-28-2015_AM'
    
    #file_names = ['I1am_Apr27_30Hz']
    #pull_filename = work_dir + 'I1_4-27-2015_AM'
    
    #file_names = ['I1pm_Apr24_30Hz']
    #pull_filename = work_dir + 'I1_4-24-2015_PM'
    
    #************ AQ2 Mouse ***********
    work_dir = '/media/cat/12TB/in_vivo/tim/yuki/AQ2/'
    
    #file_names = ['AQ2am_Jan14_30Hz']
    #pull_filename = work_dir + 'AQ2_1-14-2016_AM'
    
    file_names = ['AQ2am_Jan15_30Hz']
    pull_filename = work_dir + 'AQ2_1-15-2016_AM'
    
    lever_position_filename = work_dir + 'AQ2_1-15-2016_leverPosAM'
    
    file_names_control = []

#**********************************************************************
#************************* LOAD IMAGING DATA START ********************
#**********************************************************************

data_array=[] #Saves multiple processed data files

for file_number in range(len(file_names)):

    #Load original .tif
    images_raw = load_tif(work_dir, file_names[file_number])
    
    #Load controls for stim evoked data;
    if stimulus_evoked: images_control =  load_tif(file_names_control[file_number])
       
    #Load event times for yuki data, for example
    if False:
        trigger_value = '04'
        if yuki:  pull_times, frame_triggers, window, img_rate = load_pullfile(pull_filename, images_raw, trigger_value)

    #Load entire lever_position file:
    load_leverpositions(lever_position_filename, pull_filename, images_raw)


    #Process data
    if yuki:  
        
        if True:  #Remove 3 sec baseline
            data = remove_3sec_baseline(images_raw, frame_triggers, pull_times, work_dir, file_names[file_number])

        #Temporal filter data
        if False:
            data = filter_data()

    #Mask Ca data
    n_pixels = 128
    if True:  data = mask_data(data, work_dir, n_pixels)
    
    #Process data addtional work
    #snr_img = scipy.stats.mstats.signaltonoise(data, axis=0)
    #data = np.divide(data,snr_img)
    
    title_string = ''
    title_string = title_string+"SD Divide "
    std_img = np.std(data, axis=0)
    #data = data/std_img
    
    #print len(data)
    #print np.max(data), np.min(data)
    #data = np.log(data+np.min(data, axis=0))
    #data_ave = []
    #for k in range(len(data)):
    #    data_ave.append(np.average(data[k]))
    #data_ave =np.array(data_ave)
    #print data_ave.shape
    #data = np.divide(data, data_ave)
    
    #GENERATE ANIMATIONS
    v_max=np.ma.max(data)
    v_min=np.ma.min(data)
    #v_max = 4; v_min = -2
    print v_min, v_max

    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=5, metadata=dict(artist='Me'), bitrate=1800)

    fig = plt.figure() # make figure

    print "... generating animation..." 
    # make axesimage object
    # the vmin and vmax here are very important to get the color map correct
    im = plt.imshow(data[0], cmap=plt.get_cmap('jet'), vmin=v_min, vmax=v_max, interpolation='none')#, vmin=0, vmax=v_max)
    
    #function to update figure
    def updatefig(j):
        # set the data in the axesimage object
        im.set_array(data[j])
        plt.title(title_string +file_names[file_number] + "  #pulls: "+str(len(frame_triggers))+  
        "\nFrame: "+str(j)+"  " +str(round((float(j)-window)/img_rate,2))+"sec")
        # return the artists set
        return im,
        
    # kick off the animation
    ani = animation.FuncAnimation(fig, updatefig, frames=range(len(data)), interval=100, blit=False, repeat=True)

    if True:
    #if save_animation:
        ani.save(work_dir+file_names[file_number]+'.mp4', writer=writer)

    plt.show()

    quit()

    #********************************************************************************************
    #********************** FILTER DATA ***********************
    #********************************************************************************************
    #if (os.path.exists(work_dir+file_name+'_filtered.npy')==False):

    #*********** Use Frequency Filter ***************
    if spontaneous:
        print "Temporal filtering ..."
        from scipy.signal import butter, lfilter, filtfilt

        def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
            b, a = butter_bandpass(lowcut, highcut, fs, order=order)
            y = filtfilt(b, a, data,axis=0)
            return y

        def butter_bandpass(lowcut, highcut, fs, order=5):
            nyq = 0.5 * fs
            low = lowcut / nyq
            high = highcut / nyq
            b, a = butter(order, [low, high], btype='band')
            return b, a

        lowcut = .1
        highcut=70
        img_rate = 150        #Frame rate of imaging

        print "...band pass filtering: ", lowcut, "hz to ", highcut, "hz"
        images_raw = np.array(butter_bandpass_filter(images_raw, lowcut, highcut, img_rate, order = 2))

    if True:
        print "Spatial filtering..."
        sigma_value = 1
        images_raw = ndimage.gaussian_filter(images_raw, sigma=sigma_value) 
        
        if len(file_names_control)>0:
            images_raw_control = ndimage.gaussian_filter(images_raw_control, sigma=sigma_value)

    #*********** Use Baseline Subtraction ***************
    if True:
        if len(file_names_control)>0:
            data = (images_raw-images_raw_control)/np.average(images_raw[0:30])
        else:
            data = images_raw

    #*********** FOCUS ANALYSIS ON FRAMES 31:45:
    if 'visual' in work_dir:
        print "VISUAL clipping"
        data = data[34:84]

    #********************************************************************************************
    #********************** MASK DATA ***********************
    #********************************************************************************************
    #if (os.path.exists(work_dir+file_name+'_masked1D.npy')==False):
    if True:
        print "Masking data..."
        img_rate = 150      #Frame rate of imaging
        n_pixels = 128      #number of pixels

        #Load General mask (removes background)
        generic_mask_file = []
        generic_mask_file = work_dir + 'genericmask.txt'
        if (os.path.exists(generic_mask_file)==True):
            generic_coords = np.loadtxt(generic_mask_file)
        else:
            generic_coords = Define_generic_mask(data, work_dir, file_name, img_rate, n_pixels)

        generic_mask_indexes=np.zeros((n_pixels,n_pixels))
        for i in range(len(generic_coords)):
            generic_mask_indexes[generic_coords[i][0]][generic_coords[i][1]] = True

        temp_array = []
        #Mask all frames; NB: PROBABLY FASTER METHOD
        for i in range(0, len(data),1):
            temp = np.ma.array(data[i], mask=generic_mask_indexes, fill_value = 0., hard_mask = True)
            temp_array.append(temp)

        data = temp_array

    if True:
        print "Saving to disk..."
        np.save(work_dir + file_name+'_processed', data)

    data_array.append(data)

    #********************************************************************************************
    #********************* SHOW MOVIES *******************
    #********************************************************************************************
    if True:
        v_max=np.max(data)
        v_min=np.min(data)

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
            plt.title("Frame: "+str(j)+"\n"+str(round(float(j)/150,2))+"sec")
            # return the artists set
            return im,
        # kick off the animation
        ani = animation.FuncAnimation(fig, updatefig, frames=range(len(data)), interval=10, blit=False, repeat=True)

        if save_animation:
            ani.save(work_dir+file_name+'.mp4', writer=writer)

        plt.show()
    
    #**********************************************************
    #********************* VECTORIZE FRAMES *******************
    #**********************************************************
    #if (os.path.exists(work_dir+file_name+'_masked1D.npy')==True):

    fname_vectors = work_dir+file_name+'_masked1D'

    if True:
    #if (os.path.exists(fname_vectors+'.npy')==False):
        print "Vectorizing frames..."
        img_rate = 150      #Frame rate of imaging
        n_pixels = 128      #number of pixels

        #Load General mask (removes background)
        generic_mask_file = []
        generic_mask_file = work_dir + 'genericmask.txt'
        if (os.path.exists(generic_mask_file)==True):
            generic_coords = np.loadtxt(generic_mask_file)
        else:
            generic_coords = Define_generic_mask(data, work_dir, file_name, img_rate, n_pixels)

        generic_mask_indexes_1D=np.zeros(n_pixels*n_pixels,dtype=np.int32)
        for i in range(len(generic_coords)):
            generic_mask_indexes_1D[generic_coords[i][0]*128+ generic_coords[i][1]] = True
        
        #Flatten masked brain 2D matrix into a 1D vector
        non_masked_indexes = np.where(generic_mask_indexes_1D == 0)[0]
        b_array = []
        for i in range(len(data)):
            b_array.append(data[i].ravel()[non_masked_indexes])
        
        b_array = np.array(b_array)
        print b_array.shape
        
        fname2 = work_dir + file_name+'_dimred_'+str(n_frames)+'frames'
        print "Vectorized file: ", fname2
        if (os.path.exists(fname2+'.npy')==False):
            print "Computing corr_matrix"
            corr_matrix = np.corrcoef(b_array)
            np.save(work_dir + file_name+'_dimred_'+str(n_frames)+'frames', corr_matrix)
        else:
            corr_matrix = np.load(fname2+'.npy')
            
        if False:
            plt.imshow(corr_matrix)
            plt.show()
    
        data=corr_matrix
    else:
        print "Loading masked 1D data"
        data = np.load(fname_vectors+'.npy')
    
    #**********************************************************
    #********************* DIMENSIONALITY REDUCTION ***********
    #**********************************************************
    #Method: 0 = MDS SMACOF;  1 = t-SNE; 2 = PCA
    method = 2

    fname_dimred = work_dir + file_name+'_nodes_vertices_method_'+str(method)+'_'+str(n_frames)+'frames.npy'

    if (os.path.exists(fname_dimred)==False):
    #if True:
        ##Frame rate:
        #img_rate = 150

        ##******************** LOAD VECTORIZED DATA *********************
        #if spontaneous: 
            #length_segment = 2 * img_rate  #take first XX seconds of data
            #inner_offset = 2. * img_rate    #Offset from begining of stimulus
        #elif stimulus_evoked:
            #length_segment = len(data)
            #inner_offset = 0. * img_rate    #Offset from begining of stimulus

        #vector_data = data[inner_offset:inner_offset+length_segment]
        #print "No frames: ", len(vector_data), float(len(vector_data))/img_rate

        dim_red_data = multi_dim_scaling(data, method)

        multi_trial= []
        multi_trial.append(dim_red_data)
        single_trial=dim_red_data

        np.save(fname_dimred, multi_trial)
        
        print single_trial.shape
        print np.array(multi_trial).shape
        
    else:
        multi_trial= []
        data = np.load(fname_dimred)
        multi_trial.append(data)
        single_trial = data[0]

        print np.array(single_trial).shape
        print np.array(multi_trial).shape
        

#********************** MAKE AVERAGE VIDEO OVER ALL STIMULI **************************
if stimulus_evoked:
    import string
    fname1 = string.replace(work_dir, '/media/cat/12TB/in_vivo/tim/allen/20150303/Evoked/', '')[:-1]
    print fname1

    if (os.path.exists(work_dir+fname1+'.npy')==False):
        
        data = np.array(data_array)
        data = np.average(data_array,axis=0)
        print data.shape

        img_rate = 150      #Frame rate of imaging
        n_pixels = 128      #number of pixels

        #Load General mask (removes background)
        generic_mask_file = []
        generic_mask_file = work_dir + 'genericmask.txt'
        if (os.path.exists(generic_mask_file)==True):
            generic_coords = np.loadtxt(generic_mask_file)
        else:
            generic_coords = Define_generic_mask(data, work_dir, file_name, img_rate, n_pixels)

        generic_mask_indexes=np.zeros((n_pixels,n_pixels))
        for i in range(len(generic_coords)):
            generic_mask_indexes[generic_coords[i][0]][generic_coords[i][1]] = True

        temp_array = []
        #Mask all frames; NB: PROBABLY FASTER METHOD
        for i in range(0, len(data),1):
            temp = np.ma.array(data[i], mask=generic_mask_indexes, fill_value = 0., hard_mask = True)
            temp_array.append(temp)

        data = temp_array
        
        np.save(work_dir + fname1, data)

        #GENERATE MOVIE + SAVE
        if 'vis' in fname1:
            v_max=np.max(data[35:])
            v_min=np.min(data[35:])
            start_frame = 34
            end_frame = 104
        else:
            v_max=np.max(data)
            v_min=np.min(data)
            start_frame = 25
            end_frame = 65

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
            plt.title(fname1+"\nFrame: "+str(j)+"  "+str(round(float(j)/150,2))+"sec")
            # return the artists set
            return im,
            
        # kick off the animation
        #ani = animation.FuncAnimation(fig, updatefig, frames=range(len(data)), interval=10, blit=False, repeat=True)
        ani = animation.FuncAnimation(fig, updatefig, frames=np.arange(start_frame,end_frame,1), interval=350, blit=False, repeat=True)

        ani.save(work_dir + fname1+'.mp4', writer=writer)

        plt.show()
    
#************************* CROP/MASK TO ACTIVATED AREA ***********************
if stimulus_evoked: 
    import string
    fname1 = string.replace(work_dir, '/media/cat/12TB/in_vivo/tim/allen/20150303/Evoked/', '')[:-1]
    print fname1
    
    temp_img = np.load(work_dir + fname1+'.npy')
    print temp_img.shape
    
    #Mask it
    img_rate = 150      #Frame rate of imaging
    n_pixels = 128      #number of pixels

    #Load General mask (removes background)
    generic_mask_file = []
    generic_mask_file = work_dir + 'genericmask.txt'
    if (os.path.exists(generic_mask_file)==True):
        generic_coords = np.loadtxt(generic_mask_file)
    else:
        generic_coords = Define_generic_mask(data, work_dir, file_name, img_rate, n_pixels)

    generic_mask_indexes=np.zeros((n_pixels,n_pixels))
    for i in range(len(generic_coords)):
        generic_mask_indexes[generic_coords[i][0]][generic_coords[i][1]] = True

    temp_array = []
    #Mask all frames; NB: PROBABLY FASTER METHOD
    for i in range(0, len(temp_img),1):
        temp = np.ma.array(temp_img[i], mask=generic_mask_indexes, fill_value = 0., hard_mask = True)
        temp_array.append(temp)
    
    temp_img = temp_array#[35:]
    
    ax = plt.subplot(2,2,1)
    img_out1 = np.ma.array(np.amax(temp_img[31:51],axis=0), mask=generic_mask_indexes, fill_value = 0., hard_mask = True)
    plt.imshow(img_out1)
    plt.title("Max pixel search method")

    ax = plt.subplot(2,2,2)
    img_out2 = np.ma.array(np.average(temp_img[31:51],axis=0), mask=generic_mask_indexes, fill_value = 0., hard_mask = True)
    plt.imshow(img_out2)
    plt.title("Average pixel method")
    plt.suptitle(fname1)

    #CONTOUR MAPS
    threshold=0.4
    ax = plt.subplot(2,2,3)
    temp_img = img_out1
    temp_img = ndimage.gaussian_filter(temp_img, sigma=1)
    #temp_img = np.clip(temp_img,0,1)
    v_min = np.nanmin(temp_img)
    v_max = np.nanmax(temp_img)
    temp_img = (temp_img - v_min)/(v_max-v_min)

    temp_array = np.zeros((n_pixels,n_pixels),dtype=np.float32)
    save_indexes=[]
    for k in range(len(temp_img)):
        temp0 = np.where(np.logical_and(temp_img[k]>=threshold, temp_img[k]<=1.))[0]  #Use 1.0 for coverage maps
        save_indexes.append(temp0)
    colors=.5

    temp_array = np.zeros((n_pixels,n_pixels), dtype=np.float32)
    temp_array = np.ma.array(temp_array, mask=generic_mask_indexes, fill_value =0, hard_mask = True)
    for k in range(n_pixels):   #Scan each line
        temp_array[k][save_indexes[k]]=colors
    plt.imshow(temp_array, cmap=cm.jet)
    #plt.title("\n"+str(max_map_nspikes[i][j]),fontsize=14)

    ax = plt.subplot(2,2,4)
    temp_img = img_out2
    temp_img = ndimage.gaussian_filter(temp_img, sigma=1)
    v_min = np.nanmin(temp_img)
    v_max = np.nanmax(temp_img)
    temp_img = (temp_img - v_min)/(v_max-v_min)

    temp_array = np.zeros((n_pixels,n_pixels),dtype=np.float32)
    save_indexes=[]
    for k in range(len(temp_img)):
        temp0 = np.where(np.logical_and(temp_img[k]>=threshold, temp_img[k]<=1.))[0]  #Use 1.0 for coverage maps
        save_indexes.append(temp0)
    colors=.5

    temp_array = np.zeros((n_pixels,n_pixels), dtype=np.float32)
    temp_array = np.ma.array(temp_array, mask=generic_mask_indexes, fill_value =0, hard_mask = True)
    for k in range(n_pixels):   #Scan each line
        temp_array[k][save_indexes[k]]=colors
    plt.imshow(temp_array, cmap=cm.jet)
    #plt.title("\n"+str(max_map_nspikes[i][j]),fontsize=14)

    plt.show()

#********************** USE OPENGL TO PLOT MDS DATA **********************
print "Calling PYQT"

# prevents "The event loop is already running" messages when calling ipshell():
QtCore.pyqtRemoveInputHook()
app = QtGui.QApplication(sys.argv) #Cat: app start

mainwindow = MainWindow()

print "Done INITIALIZATION"
mainwindow.work_dir = work_dir
mainwindow.file_names = file_names
mainwindow.single_trial = single_trial
mainwindow.multi_trial = multi_trial
mainwindow.time_segs = len(single_trial)
mainwindow.method = method

mainwindow.show()

sys.exit(app.exec_())


#*********************************************************************

#Matplotlib 3d 
#fig = plt.figure()
#ax = fig.gca(projection='3d')

#multi_trials=[]

##Compute MDS reduction individually for each movie repeat
#if False:
    #for i in range(3):
    ##for i in range(len(cell_rate_trans)):

        #print "Movie repeat: ", i
        #dists = sklearn.metrics.pairwise.pairwise_distances(cell_rate_trans[i])

        #adist = np.array(dists)
        #amax = np.amax(adist)
        #adist /= amax

        #mds = manifold.MDS(n_components=3, dissimilarity="euclidean", random_state=6)
        #results = mds.fit(adist)

        #coords = results.embedding_

        #if plotting:
            #for i in range(len(coords)-1):
                #loc_1 = [coords[i,0],coords[i,1],coords[i,2]]
                #loc_2 = [coords[i+1,0],coords[i+1,1],coords[i+1,2]]

                #ax.plot([loc_1[0],loc_2[0]],[loc_1[1],loc_2[1]],[loc_1[2],loc_2[2]],color=(float(i)/len(coords),0,float(i)/len(coords)))
                #ax.scatter(coords[i, 0],coords[i, 1],coords[i, 2], s=20, color=(float(i)/len(coords),0,float(i)/len(coords)))
        #multi_trials.append(coords*1E3)

    #if plotting:
        #plt.show()













