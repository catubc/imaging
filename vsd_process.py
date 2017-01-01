from tsf_ptcs_classes import *
from sequential_firing import *
from opengl_pyqt_classes import *

import numpy as np
import matplotlib.pyplot as plt
import struct
import glob
from scipy import ndimage
import matplotlib.animation as animation

from PIL import Image


#**********************************
#work_dir = '/media/cat/12TB/in_vivo/tim/alex/2016-03-10_FSF/Active Stimuli/'
#work_dir = '/media/cat/12TB/in_vivo/tim/alex/2016-03-10_FSF/'

#work_dir = '/media/cat/12TB/in_vivo/tim/alex/2016-03-14/FSF/'
#work_dir = '/media/cat/12TB/in_vivo/tim/alex/2016-03-14/led/'


#work_dir = '/media/cat/12TB/in_vivo/tim/alex/2016-03-15/FSF_0iso/'
#work_dir = '/media/cat/12TB/in_vivo/tim/alex/2016-03-15/FSF_1iso/'
#work_dir = '/media/cat/12TB/in_vivo/tim/alex/2016-03-15/FSF_1iso_later/'

#tifs = glob.glob(work_dir+ "*hi*.tif")

#****************** SPONTANEOUS ACTIVITY *****************
work_dir = '/media/cat/12TB/in_vivo/tim/alex/2016-03-15/spont_1iso/'
#work_dir = '/media/cat/12TB/in_vivo/tim/alex/2016-03-14/spont/'
tifs = glob.glob(work_dir+ "*.tif")

data_array=[]
data_array_control=[]
data_out_array=[]
#Load .tif files
for i in [0]:#range(len(tifs)):
    
    #if i >5: continue
    
    if 'spont' in work_dir:
        fname = tifs[i][:-4]
    else:
        fname = work_dir+'hi'+str(i+1)
    if (os.path.exists(fname+'.npy')==False):
        img = Image.open(fname+'.tif')

        counter=0
        if True:
            while True:
                try:
                    img.seek(counter)
                except EOFError:
                    break
                counter+=1
                #print counter

        n_pixels = 128
        images_raw = np.zeros((counter, n_pixels, n_pixels), dtype = np.float32)

        #print "n_frames: ", n_frames
        for p in range(0, counter,1): 
            try:
                img.seek(p)
                images_raw [p] = img 

            except EOFError:
                break

        print "Saving imaging array..."

        np.save(fname, images_raw)

    else:
        print "Loading from .npy"
        images_raw = np.load(fname+'.npy')

    #Load control
    if 'hi' in tifs[i]:
        fname = work_dir+'no'+str(i+1)
        if (os.path.exists(fname+'.npy')==False):
            img = Image.open(fname+'.tif')

            counter=0
            if True:
                while True:
                    try:
                        img.seek(counter)
                    except EOFError:
                        break
                    counter+=1

            n_pixels = 128
            images_raw_control = np.zeros((counter, n_pixels, n_pixels), dtype = np.float32)

            print "n_frames: ", counter
            for p in range(0, counter,1): 
                try:
                    img.seek(p)
                    images_raw_control[p] = img 

                except EOFError:
                    break                    
        
            np.save(fname, images_raw_control)

        else:
            
            images_raw_control = np.load(fname+'.npy')
    
    #Remove baseline:
    #images_raw = ndimage.gaussian_filter(images_raw, sigma=1) 
    #images_raw_control = ndimage.gaussian_filter(images_raw_control, sigma=1) 

    #baseline = np.average(images_raw[0:30])
    #temp_data = []
    #for k in range(len(images_raw)):
    #    temp_data.append((images_raw[k] - images_raw_control[k])/np.average(images_raw[0:30], axis=0))
    #data = ndimage.gaussian_filter(data, sigma=1) 
    
    baseline = np.average(images_raw, axis=0)
    data_out = images_raw/baseline
    data_out_array.append(data_out)

    if 'hi' in tifs[i]:
        data_array.append(images_raw)
        data_array_control.append(images_raw_control)

    if 'spont' in work_dir: break


#if 'hi' in tifs[i]:
    #data = np.average(data_array, axis=0)
    #data_control = np.average(data_array_control, axis=0)

    #data = (data-data_control)/np.average(data[0:30])
    #data = ndimage.gaussian_filter(data, sigma=1) 

#data = np.average(data_out_array, axis=0)
#data = ndimage.gaussian_filter(data_out, sigma=1) 

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
#data_array = np.average(data_array, axis=0)

data = np.array(butter_bandpass_filter(data_out[:], lowcut, highcut, img_rate, order = 2))
data = ndimage.gaussian_filter(data, sigma=0.75)[0:2000]


if True:
    #v_max=np.max(data)
    #v_min=np.min(data)

    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=5, metadata=dict(artist='Me'), bitrate=1800)

    fig = plt.figure() # make figure

    # make axesimage object
    # the vmin and vmax here are very important to get the color map correct
    im = plt.imshow(data[0], cmap=plt.get_cmap('jet'), interpolation='none')#, vmin=0, vmax=v_max)
    # function to update figure
    def updatefig(j):
        # set the data in the axesimage object
        im.set_array(data[j])
        plt.title("Frame: "+str(j)+"\n"+str(round(float(j)/150,2))+"sec")
        # return the artists set
        return im,
    # kick off the animation
    ani = animation.FuncAnimation(fig, updatefig, frames=range(len(data)), interval=10, blit=False, repeat=True)

    if False:
        ani.save(work_dir+'vid.mp4', writer=writer)

    plt.show()
        
        

quit()
