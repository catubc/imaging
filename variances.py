import matplotlib.pyplot as plt
import numpy as np
import os, csv
from sklearn import cluster, datasets
import glob
import scipy


def quick_mask_allframe(data, midline_mask_n_pixels):
    
    n_pixels = len(data[0])
    print "...n_pixels: ", n_pixels
    
    generic_mask_file = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-11/genericmask.txt'
    if (os.path.exists(generic_mask_file)==True):
        generic_coords = np.int32(np.loadtxt(generic_mask_file))
    else:
        print "...generic mask not found..."
        return
    
    #Load generic mask
    generic_mask_indexes=np.zeros((256,256))
    for i in range(len(generic_coords)):
        generic_mask_indexes[int(generic_coords[i][0])][int(generic_coords[i][1])] = True

    if n_pixels != 256: #Subsample mask:
        generic_mask_indexes = scipy.misc.imresize(generic_mask_indexes,n_pixels/256.)

    #Load midline mask
    for i in range(midline_mask_n_pixels):
        generic_mask_indexes[:,n_pixels/2+int(midline_mask_n_pixels/2)-i]=True

    #Apply full mask; probably FASTER METHOD
    temp_array = np.ma.array(np.zeros((len(data), n_pixels,n_pixels),dtype=np.float32), mask=True)
    for k in range(len(data)):
        temp_array[k] = np.ma.masked_array(data[k], mask=generic_mask_indexes, fill_value=0)
    
    return temp_array


#************************************************************
if True: 
    #file_name= '/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-2/img_avg_2015-7-22-2_unit00_ch05_all_3sec_window_00965_spikes.npy'
    #file_name = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-8-9electrodein-iso0/img_avg_2015-11-18-8-9electrodein-iso0_unit01_ch11_all_3sec_window_01476_spikes.npy'
    #file_name = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/img_avg_2015-11-27-16-deep-iso0_unit14_ch12_all_3sec_window_03360_spikes.npy'
    file_name = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-15-allelectrodeinthalamus-is0/img_avg_2015-12-2-15-allelectrodeinthalamus-is0_unit21_ch09_all_3sec_window_01445_spikes.npy'
    
    data = np.load(file_name)

    data = quick_mask_allframe(data, 3)

    v_max = np.nanmax(np.abs(data))
    print v_max

    for k in range(len(data)):
        print "...map: ", k
        ax = plt.subplot(10, 18, k+1)
        #ax = plt.subplot(16, 20, k+1)
        plt.imshow(data[k], vmin=-v_max, vmax = v_max)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

    plt.suptitle(file_name+" DF/F: "+str(round(v_max*100,1))+"%", fontsize=15)
    plt.show()


#************************************************************
#file_name = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-2/stack1D_2015-7-22-2_unit00_ch05_all_3sec_window_00965_spikes.npy'
#file_name = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-8-9electrodein-iso0/stack1D_2015-11-18-8-9electrodein-iso0_unit01_ch11_all_3sec_window_01476_spikes.npy'
#file_name = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/stack1D_2015-11-27-16-deep-iso0_unit14_ch12_all_3sec_window_03360_spikes.npy'
file_name = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-15-allelectrodeinthalamus-is0/stack1D_2015-12-2-15-allelectrodeinthalamus-is0_unit21_ch09_all_3sec_window_01445_spikes.npy'

print "...loading 1D stack (large file)..."
data = np.load(file_name)

print "...computing variances..."
data = np.var(data, axis=0)
    
img_stack = []
for p in range(data.shape[1]/64):
    img_stack.append(data[:,p*64:(p+1)*64])

data= np.array(img_stack)

data = quick_mask_allframe(data, 3)

v_max = np.nanmax(np.abs(data))
print v_max

for k in range(len(data)):
    print "...map: ", k
    ax = plt.subplot(10, 18, k+1)
    #ax = plt.subplot(16, 20, k+1)
    plt.imshow(data[k], vmin=-v_max, vmax = v_max)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

plt.suptitle(file_name+" DF/F: "+str(round(v_max*100,1))+"%", fontsize=15)
plt.show()





