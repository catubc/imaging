import os
import glob
import numpy as np
import struct
import string, re
import scipy
import tifffile as tiff
from skimage.measure import block_reduce
import matplotlib.pyplot as plt
import shutil
import matplotlib.patches as mpatches

from scipy import ndimage
    
def quick_mask_single_frame(data, midline_mask_n_pixels):
    
    n_pixels = 256
        
    generic_mask_file = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-11/genericmask.txt'
    if (os.path.exists(generic_mask_file)==True):
        generic_coords = np.int32(np.loadtxt(generic_mask_file))
    else:
        print "...generic mask not found..."
        return
    
    #Load generic mask
    generic_mask_indexes=np.zeros((n_pixels,n_pixels))
    for i in range(len(generic_coords)):
        generic_mask_indexes[int(generic_coords[i][0])][int(generic_coords[i][1])] = True

    #Load midline mask
    for i in range(midline_mask_n_pixels):
        generic_mask_indexes[:,n_pixels/2+int(midline_mask_n_pixels/2)-i]=True

    #Apply full mask; probably FASTER METHOD
    temp_array = np.ma.array(np.zeros((n_pixels,n_pixels),dtype=np.float32), mask=True)
    temp_array = np.ma.masked_array(data, mask=generic_mask_indexes, fill=0)
    
    return temp_array
    
    
#****************************************************************************************************


#data = np.load('/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-11/2015-12-11-5-10electincortex-iso0_lfp/2015-12-11-5-10electincortex-iso0_lfp_lfpmap_band_alpha_channel_15.npy')
#plt.imshow(data)
#plt.show()
#quit()

bands = ['delta', 'theta', 'alpha', 'beta', 'gamma']#, 'high']


#base_name = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-11/2015-12-11-5-10electincortex-iso0_lfp/2015-12-11-5-10electincortex-iso0_lfp_lfpmap_band' #_delta_channel_00.npy'
#base_name = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-11/2015-12-11-5-10electincortex-iso0_lfp/2015-12-11-5-10electincortex-iso0_lfp_powermap_band' #_delta_channel_00.npy'

#base_name = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-11/2015-12-11-3-10electincortex-iso1_lfp/2015-12-11-3-10electincortex-iso1_lfp_lfpmap_band'
#base_name = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-11/2015-12-11-16-allelectinthalamus-iso0_lfp/2015-12-11-16-allelectinthalamus-iso0_lfp_lfpmap_band'

base_name = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-14-allelectrodeinthalamus-is0_lfp/2015-12-2-14-allelectrodeinthalamus-is0_lfp_lfpmap_band'



mid_mask = [5,10,15,25,10,20]
max_percentiles = [99.5,99.5,99.5,99.5,99.5,97]

fig, axes = plt.subplots(nrows=4, ncols=4)
fig.tight_layout()

chs = np.arange(0,16,1)

sigma = 3

for ctr, band in enumerate(bands):
    print band
    for ch in chs:
        #ax = plt.subplot(16,6, (ch-6)*6+ctr+1)
        ax = plt.subplot(16,6, ch*6+ctr+1)
        data = np.load(base_name+'_'+band+'_channel_'+str(ch).zfill(2)+'.npy')
        data = ndimage.gaussian_filter(data, sigma=sigma)
                
        data = quick_mask_single_frame(data, mid_mask[ctr])

        #Make single trial motion mask:
        

        #v_max = np.nanmax(np.abs(data)); v_min = -v_max
        #v_max = np.nanmax(data); v_min = np.nanmin(data)
        max_percentile = max_percentiles[ctr]
        min_percentile = 1
        data_1d = data.ravel()
        data_1d = data_1d[np.logical_not(np.isnan(data_1d))]        #Remove 
        v_max = abs(np.percentile(data_1d, max_percentile))               #Mark stroke as the 97.5 percentile and higher values; 
        v_min = np.percentile(data_1d, min_percentile)               #Mark stroke as the 97.5 percentile and higher values; 
        #print "...masked percentile value: ", v_max
        #v_max = v_max; v_min = -v_max
        print v_min, v_max
        


        plt.imshow(data, vmin=v_min, vmax=v_max)
        #plt.title(band + "  " + str(ch))
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

plt.show()

plt.suptitle(file_name)
plt.ylim(-1, 1)
plt.show()


quit()
sum_temp = []
for k in range(len(data)):
    line = np.zeros(128, dtype=np.int32)
    line[data[k]]=1.0
    sum_temp.append(line)
        
sum_temp = np.vstack((sum_temp))
plt.imshow(sum_temp)
plt.show()
