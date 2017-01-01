#from tsf_ptcs_classes import *
#from sequential_firing import *
#from opengl_pyqt_classes import *
#from sta_utils import *

from mouse import *

from lever_pull_utils import *
import filter
import inspect
import os
import glob
import multiprocessing as mp

import numpy as np
import matplotlib.pyplot as plt
import struct
import string
import scipy
import skimage
from skimage import data
from skimage.transform import rotate

from matplotlib import cm
import matplotlib.animation as animation
from scipy import ndimage


#Initialize pyqtgraph - needs to be done at beginning - or even better in a main conditional below...
from pyqtgraph.Qt import QtCore, QtGui

app = QtGui.QApplication([])    #Starts up the QtGui; makes gl object that needs to be passed to graph routines...


#**********************************
#************* DATA FILES *********
#**********************************
home_dir = '/media/cat/12TB/in_vivo/tim/yuki/'

#GCamp6f mice; NOT good overall; odd responses; 
#mice_news=['AQ2', 'AQ3', 'AQ5', 'AV2', 'I1', 'T1', 'T2']

#GCamp6s mice
#mice_names =['IA1', 'IA2', 'IA3', 'IJ1', 'IJ2']
#mice_names = ['AI3','AK4','AK5','AR4','BA1','BA2','Y1']
mice_names = ['IA2']

 
#**********************************
#********** LOAD MICE  ************
#**********************************
n_sec = 3   #number of sec for window to compute DF/F

mice = []   #list of mice to load
for mouse_name in mice_names:
    if False:
        mouse = Mouse(mouse_name, home_dir, n_sec)#, window, n_pixels)
        mouse.process()
    else:
        mouse = Mouse(mouse_name, home_dir, n_sec)#, window, n_pixels)
        mouse.load()


#*************************************
#*** PARSE DATA INTO CHUNK PERIODS ***
#*************************************
pre_stroke_days = 21;  post_stroke_days = 14;  post_post_stroke_days = 42
mouse.chunk_data(pre_stroke_days, post_stroke_days, post_post_stroke_days)


#*************************************
#****** DIM REDUCTION OF TRACES ******
#*************************************
#methods: 0: MDS; 1:tSNE; 2: PCA; 3: Barnes-Hut tSNE
dim_red_data = multi_dim_scaling(mouse, 2) 


#*************************************
#************* CLUSTERING ************
#*************************************
#Clustering: Kmeans vs. Meanshift
mouse.n_clusters = 81
mouse.cluster_labels = KMEANS(dim_red_data, mouse.n_clusters)

#Plot 3d distributions 
if False: plot_pyqt(app, dim_red_data, mouse.cluster_labels)

#Plot 2d distributions
plot_traces_pyqt(app, mouse)


counter=0
for p in range(len(mouse.traces)): 
    chunk_labels = mouse.cluster_labels[counter:counter+len(mouse.DFF_files[p])]
    chunk_indexes =  np.where(chunk_labels==mouse.selected_cluster)[0]
    counter+=len(mouse.DFF_files[p])

    #Select traces and DFF data for selected cluster only
    mouse.traces[p] = mouse.traces[p][chunk_indexes]
    mouse.DFF_files[p] = mouse.DFF_files[p][chunk_indexes]
    mouse.trial_times[p] = mouse.trial_times[p][chunk_indexes]
    mouse.days_counter[p] = mouse.days_counter[p][chunk_indexes]

    print len(mouse.traces[p])


#*************************************
#******** PLOT [Ca] DYNAMICS *********
#*************************************

#PLOT AVERAGE [Ca] DYNAMICS  
#plot_selected_cluster_DFF(mouse)

plot_single_trial(mouse)


print "Clean exit..."

quit()




#LOAD MASK
generic_mask_indexes = load_generic_mask(main_dir, data_array)

#CHECK DATA FOR BAD LEVER PULLS 
if False: 
    bad_trials=[]
    for p in range(len(data_array)):
        print p
        if np.max(data_array[p])>.5: 
            print "bad trial"
            bad_trials.append(p)
            
    print bad_trials
    data_array_good_trials = np.delete(data_array, bad_trials, axis=0)
    print data_array_good_trials.shape
    np.save(home_dir+mouse+"/"+mouse+"_imaging_good.npy", data_array_good_trials)

    all_traces_good_trials = np.delete(all_traces, bad_trials, axis=0)
    all_traces = np.save(home_dir+mouse+"/"+mouse+"_traces_good.npy", all_traces_good_trials)


#PLOT INDIVIDUAL TRIALS
start_trial = 800
end_trial = 810
trials = np.arange(start_trial, end_trial, 1)

if True:
    #Plot all data absolute scale
    data_plot = []
    for p in trials:
        masked = np.ma.array(np.zeros((120,n_pixels,n_pixels),dtype=np.float32), mask=True)
        masked = fast_mask(data_array[p], generic_mask_indexes)
        rows=[]
        for k in range(20):
            ave_plot = np.ma.average(masked[k*6:(k+1)*6], axis=0)
            rows.append(ave_plot)
        data_plot.append(np.ma.hstack((rows)))
    data_plot = np.ma.vstack((data_plot))
    
    v_min = np.min(data_plot)
    v_max = np.max(data_plot)
    
    ax1=plt.subplot(211)
    plt.imshow(data_plot, vmin=v_min, vmax=v_max)
    ax1.yaxis.set_ticks([])
    ax1.xaxis.set_ticks([])
    plt.ylabel("Trials", fontsize=20)
    plt.title("Absolute normalized, DF/F min:"+str(int(v_min*100))+ "%, max: "+str(int(v_max*100)) + "%",  fontsize=20 )
    #plt.show()
    
    #Plot all data relative scale
    data_plot = []
    for p in trials:
        masked = np.ma.array(np.zeros((120,n_pixels,n_pixels),dtype=np.float32), mask=True)
        masked = fast_mask(data_array[p], generic_mask_indexes)
        rows=[]
        for k in range(20):
            ave_data = np.ma.average(masked[k*6:(k+1)*6], axis=0)
            ave_plot = (ave_data-np.min(ave_data))/(np.max(ave_data)-np.min(ave_data))
            rows.append(ave_plot)
        data_plot.append(np.ma.hstack((rows)))
    data_plot = np.ma.vstack((data_plot))
        
    ax1=plt.subplot(212)
    plt.imshow(data_plot)#, vmin=v_min, vmax=v_max)
    ax1.yaxis.set_ticks([])
    #ax1.xaxis.set_ticks([])
    print data_plot.shape
    
    old_xlabel = np.linspace(0,data_plot.shape[1], 20)
    new_xlabel = np.around(np.linspace(-2.0, 2.0, 20), decimals=1)
    plt.xticks(old_xlabel, new_xlabel)
    
    plt.title("Frame normalized", fontsize=20)
    plt.xlabel("Time from reward code (sec)", fontsize = 20)
    plt.suptitle(mouse + " trials: "+str(start_trial)+ " .. "+ str(end_trial), fontsize = 30)
    plt.ylabel("Trials", fontsize=20)

    plt.show()

quit()

#PLOT AVERAGES
n_procs = 24
n_trials = 5 #len(data_array) #Number of trials in each row of averates
#data_array = data_array[0:2000]
#for p in range(0, len(data_array), n_trials):
for p in range(1040, 1080, n_trials):
    print "...plotting trials section: ", p
    print "...averaging data..."
    
    
    data_plot = np.average(data_array[p:p + n_trials],axis=0)
    #data_plot = np.average(data_array[p:n_trials-4000],axis=0)
    
    
    #data_chunks = []
    #for i in range(0,len(data_array), len(data_array)/n_procs):
        #data_chunks.append(data_array[i:i+len(data_array)/n_procs])
    
    #data_chunks=np.array(data_chunks[:-1])
    #print data_chunks.shape
    
    #data_plot=[]
    #pool = mp.Pool(n_procs+1)
    #data_plot.extend(pool.map(average_parallel,data_chunks))
    
    #pool.close()
   
    #data_plot = np.average(data_plot, axis=0) #This should be ok as the chunks are roughly the same size...
    
    #plt.imshow(data_plot[0])
    #plt.show()
    
    data_plot = mask_data(data_plot, main_dir)
    
    print "...plotting data..."
    v_min=0
    v_max=0
    for k in range(20):
        ax=plt.subplot(10,20,k+1)
        chunk_ave = np.ma.average(data_plot[k*6:(k+1)*6], axis=0) 
        plt.imshow(chunk_ave )#, vmin=-0.05, vmax=0.05)
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        
        #if k==0: ax.set_ylabel("#: "+str(p+1))
        if (p%n_trials)==0: ax.set_title(str(float(k)/5-2.0)+"s")
        v_min = min(v_min, np.min(chunk_ave))
        v_max = max(v_max, np.max(chunk_ave))
  

    print "...plotting data..."
    print v_min, v_max
    for k in range(20):
        ax=plt.subplot(10,20,20+k+1)
        plt.imshow(np.ma.average(data_plot[k*6:(k+1)*6], axis=0), vmin=v_min, vmax=v_max)
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        
        #if k==0: ax.set_ylabel("#: "+str(p+1))
        #if (p%n_trials)==0: ax.set_title(str(float(k)/5-2.0)+"s")

        if k==0: plt.text(-70, -140, str(p)+ " .. " + str(p+n_trials), rotation=90, fontsize=20)
    plt.suptitle(mouse + " trials: "+ str(p)+ " .. " + str(p+n_trials), fontsize=25)
    plt.show()



quit()
quit()




#************** COMBINE ALL DATA *******************
if False:
    pass

#************** SHOW INDIVIDUAL TRIALS ************
else:
    print "GENERATING VIDS..."
    #print len(good_tifs), len(data_array_individual)
    
    vid_array = []
    n_pulls = []
    tif_names = []
    for k in range(len(data_array_individual)):
        print "...masking vid: ", k
        if len(data_array_individual[k])==0: continue   #if there are zero correct pulls skip
        
        temp_array = np.array(data_array_individual[k])
        if temp_array.shape[0]<10: continue #Exit if less than 10 pulls
        n_pulls.append(temp_array.shape[0])

        temp_vid = np.average(temp_array, axis=0)    #Average over all trials recorded from

        n_pixels = 128
        if True:  temp_vid = mask_data(temp_vid, home_dir + mouse, n_pixels)

        vid_array.append(temp_vid)
        
        tif_names.append(good_tifs[k])
        
        if len(vid_array)==40: break

    #***********GENERATE ANIMATIONS
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=5, metadata=dict(artist='Me'), bitrate=1800)

    fig = plt.figure() # make figure

    print "... generating animation..." 
    # make axesimage object
    # the vmin and vmax here are very important to get the color map correct

    im=[]
    for k in range(len(vid_array)):
        ax = plt.subplot(5,8,k+1)
        
        #v_max=np.ma.max(vid_array[k])
        #v_min=np.ma.min(vid_array[k])
        v_max=0.01
        v_min=-0.01
        
        #print v_max, v_min
        
        ax.get_xaxis().set_visible(False)
        ax.yaxis.set_ticks([])
        #ax.set_axislabel('Galactic Longitude', minpad=0.3)
        ax.yaxis.labelpad = 0
        ax.set_ylabel(tif_names[k][:-9],fontsize=6)
        ax.set_title("#: "+str(n_pulls[k]))
        
        im.append([])
        im[k] = plt.imshow(vid_array[k][0], cmap=plt.get_cmap('jet'), vmin=v_min, vmax=v_max, interpolation='none')#, vmin=0, vmax=v_max)

    #function to update figure
    def updatefig(j):
        print "... updating frame: ", j
        plt.suptitle(mouse + ",  reward code: "+ reward_code+ ", v_max: "+str(int(v_max*100))+ "%  v_min: "+str(int(v_min*100))+"%"+
        "\nFrame: "+str(j)+"  " +str(round((float(j)-window)/img_rate,2))+"sec")

        # set the data in the axesimage object
        for k in range(len(vid_array)):
            im[k].set_array(vid_array[k][j])

        # return the artists set
        return im
        
    # kick off the animation
    ani = animation.FuncAnimation(fig, updatefig, frames=range(len(vid_array[0])), interval=100, blit=False, repeat=True)

    if True:
    #if save_animation:
        ani.save(home_dir+mouse+"/"+mouse+"_reward_code_"+reward_code+'.mp4', writer=writer)

    plt.show()



#************* PROCESS DATA **************
#Process data addtional work
#import scipy
#snr_img = scipy.stats.mstats.signaltonoise(data, axis=0)
#data = np.divide(data,snr_img)

title_string = ''
#title_string = title_string+"SD Divide "
#std_img = np.std(data, axis=0)
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



quit()

if True:
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






#NORMALITY TEST
def normality_test():
    ax=plt.subplot(211)
    for k in range(0, len(all_traces), n_trials):
        trace_ave = np.average(all_traces[k:k +n_trials], axis=0)
        plt.plot(t, trace_ave, c=cm.jet(int(k/n_trials*color_scale)), linewidth=3, alpha=1)
        temp_std = np.std(all_traces[k:k +n_trials], axis=0)
        #if k == 0: ax.fill_between(t, trace_ave+temp_std, trace_ave-temp_std, color=cm.jet(int(k/n_traces*color_scale)), alpha=0.3)

    #ax.fill_between(t, trace_ave+temp_std, trace_ave-temp_std, color=cm.jet(int(k/n_trials*color_scale)), alpha=0.3)

    reward_code = '04'
    plt.suptitle(mouse+ ",  code: "+ reward_code + ",  # pulls: "+str(len(all_traces)), fontsize =30)
    plt.ylabel("Lever Angle (and std)", fontsize=25)
    #plt.xlabel("Time from code threshold (sec)", fontsize=25)
    plt.ylim(-5,85)
    plt.xlim(t[0],t[-1])
    plt.plot([0,0],[-25,85], 'r--', linewidth=2, color='black', alpha=0.4)
    plt.plot([t[0],t[-1]],[0,0], 'r--', linewidth=2, color='black', alpha=0.4)
    #plt.show()
    
    
    ax=plt.subplot(212)
    p_val_norm = []
    for k in range(100):
        p_val_norm.append(scipy.stats.mstats.normaltest(all_traces[:,k])[1]*100., axis=0)

    plt.plot(t, p_val_norm, linewidth=2)
    plt.ylabel("Probability that Gaussian", fontsize=25)
    plt.xlabel("Time from code threshold (sec)", fontsize=25)
    plt.ylim(-1,100)
    plt.xlim(t[0],t[-1])
    plt.plot([0,0],[0,100], 'r--', linewidth=2, color='black', alpha=0.4)
    plt.plot([t[0],t[-1]],[0,0], 'r--', linewidth=2, color='black', alpha=0.4)

    mu, sigma = 0, .01 # mean and standard deviation
    s = np.random.normal(mu, sigma, 10000)
    print scipy.stats.mstats.normaltest(s, axis=0)
    
    plt.show()
    
    quit()







