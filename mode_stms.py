import matplotlib.pyplot as plt
import numpy as np
import os, csv
from sklearn import cluster, datasets
import glob
import scipy


def quick_mask_single_allframe(data, midline_mask_n_pixels):
    
    n_pixels = len(data[0])
    #print "...n_pixels: ", n_pixels
    
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
        generic_mask_indexes = scipy.misc.imresize(generic_mask_indexes,.25)

    #Load midline mask
    for i in range(midline_mask_n_pixels):
        generic_mask_indexes[:,n_pixels/2+int(midline_mask_n_pixels/2)-i]=True

    #Apply full mask; probably FASTER METHOD
    temp_array = np.ma.array(np.zeros((len(data), n_pixels,n_pixels),dtype=np.float32), mask=True)
    for k in range(len(data)):
        temp_array[k] = np.ma.masked_array(data[k], mask=generic_mask_indexes, fill_value=0)
    
    return temp_array


#******************************************************************************

files=[
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-8-9electrodein-iso0/img_avg_2015-11-18-8-9electrodein-iso0_unit03_ch13_modes_3sec_window_02076_spikes_first.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-8-9electrodein-iso0/img_avg_2015-11-18-8-9electrodein-iso0_unit03_ch13_modes_3sec_window_02076_spikes_tonic.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-8-9electrodein-iso0/img_avg_2015-11-18-8-9electrodein-iso0_unit03_ch13_modes_3sec_window_02076_spikes_last.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-8-9electrodein-iso0/img_avg_2015-11-18-8-9electrodein-iso0_unit03_ch13_modes_3sec_window_02076_spikes_burst.npy'

#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-11/2015-12-11-12-allelectinthalamus-iso1/img_avg_2015-12-11-12-allelectinthalamus-iso1_unit02_ch03_modes_3sec_window_01213_spikes_first.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-11/2015-12-11-12-allelectinthalamus-iso1/img_avg_2015-12-11-12-allelectinthalamus-iso1_unit02_ch03_modes_3sec_window_01213_spikes_tonic.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-11/2015-12-11-12-allelectinthalamus-iso1/img_avg_2015-12-11-12-allelectinthalamus-iso1_unit02_ch03_modes_3sec_window_01213_spikes_last.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-11/2015-12-11-12-allelectinthalamus-iso1/img_avg_2015-12-11-12-allelectinthalamus-iso1_unit02_ch03_modes_3sec_window_01213_spikes_burst.npy'

#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/img_avg_2015-11-27-16-deep-iso0_unit01_ch07_modes_3sec_window_04807_spikes_first.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/img_avg_2015-11-27-16-deep-iso0_unit01_ch07_modes_3sec_window_04807_spikes_tonic.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/img_avg_2015-11-27-16-deep-iso0_unit01_ch07_modes_3sec_window_04807_spikes_last.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/img_avg_2015-11-27-16-deep-iso0_unit01_ch07_modes_3sec_window_04807_spikes_burst.npy'

#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/img_avg_2015-11-27-16-deep-iso0_unit00_ch06_modes_3sec_window_02552_spikes_first.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/img_avg_2015-11-27-16-deep-iso0_unit00_ch06_modes_3sec_window_02552_spikes_tonic.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/img_avg_2015-11-27-16-deep-iso0_unit00_ch06_modes_3sec_window_02552_spikes_last.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/img_avg_2015-11-27-16-deep-iso0_unit00_ch06_modes_3sec_window_02552_spikes_burst.npy',

#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/img_avg_2015-11-27-16-deep-iso0_unit08_ch11_modes_3sec_window_01419_spikes_first.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/img_avg_2015-11-27-16-deep-iso0_unit08_ch11_modes_3sec_window_01419_spikes_tonic.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/img_avg_2015-11-27-16-deep-iso0_unit08_ch11_modes_3sec_window_01419_spikes_last.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/img_avg_2015-11-27-16-deep-iso0_unit08_ch11_modes_3sec_window_01419_spikes_burst.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/img_avg_2015-11-27-16-deep-iso0_unit08_ch11_all_3sec_window_01419_spikes.npy'

#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/img_avg_2015-11-27-16-deep-iso0_unit14_ch12_modes_3sec_window_03360_spikes_first.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/img_avg_2015-11-27-16-deep-iso0_unit14_ch12_modes_3sec_window_03360_spikes_tonic.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/img_avg_2015-11-27-16-deep-iso0_unit14_ch12_modes_3sec_window_03360_spikes_last.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/img_avg_2015-11-27-16-deep-iso0_unit14_ch12_modes_3sec_window_03360_spikes_burst.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/img_avg_2015-11-27-16-deep-iso0_unit14_ch12_all_3sec_window_03360_spikes.npy'

#OK CELL, BUT modes do not cluster;
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-13-deep-iso0/img_avg_2015-11-18-13-deep-iso0_unit14_ch10_modes_3sec_window_02140_spikes_first.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-13-deep-iso0/img_avg_2015-11-18-13-deep-iso0_unit14_ch10_modes_3sec_window_02140_spikes_tonic.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-13-deep-iso0/img_avg_2015-11-18-13-deep-iso0_unit14_ch10_modes_3sec_window_02140_spikes_last.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-13-deep-iso0/img_avg_2015-11-18-13-deep-iso0_unit14_ch10_modes_3sec_window_02140_spikes_burst.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-13-deep-iso0/img_avg_2015-11-18-13-deep-iso0_unit14_ch10_all_3sec_window_02140_spikes.npy'

#OK Burst vs All spikes
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-14-deep-iso0/img_avg_2015-11-18-14-deep-iso0_unit07_ch06_modes_3sec_window_06810_spikes_first.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-14-deep-iso0/img_avg_2015-11-18-14-deep-iso0_unit07_ch06_modes_3sec_window_06810_spikes_tonic.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-14-deep-iso0/img_avg_2015-11-18-14-deep-iso0_unit07_ch06_modes_3sec_window_06810_spikes_last.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-14-deep-iso0/img_avg_2015-11-18-14-deep-iso0_unit07_ch06_modes_3sec_window_06810_spikes_burst.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-14-deep-iso0/img_avg_2015-11-18-14-deep-iso0_unit07_ch06_all_3sec_window_06810_spikes.npy'

#OK all vs. first; burst shifted negatively in time
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-14-deep-iso0/img_avg_2015-11-18-14-deep-iso0_unit35_ch14_modes_3sec_window_01052_spikes_first.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-14-deep-iso0/img_avg_2015-11-18-14-deep-iso0_unit35_ch14_modes_3sec_window_01052_spikes_tonic.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-14-deep-iso0/img_avg_2015-11-18-14-deep-iso0_unit35_ch14_modes_3sec_window_01052_spikes_last.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-14-deep-iso0/img_avg_2015-11-18-14-deep-iso0_unit35_ch14_modes_3sec_window_01052_spikes_burst.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-14-deep-iso0/img_avg_2015-11-18-14-deep-iso0_unit35_ch14_all_3sec_window_01052_spikes.npy'

#OK all vs first;
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-14-deep-iso0/img_avg_2015-11-18-14-deep-iso0_unit38_ch15_modes_3sec_window_03893_spikes_first.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-14-deep-iso0/img_avg_2015-11-18-14-deep-iso0_unit38_ch15_modes_3sec_window_03893_spikes_tonic.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-14-deep-iso0/img_avg_2015-11-18-14-deep-iso0_unit38_ch15_modes_3sec_window_03893_spikes_last.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-14-deep-iso0/img_avg_2015-11-18-14-deep-iso0_unit38_ch15_modes_3sec_window_03893_spikes_burst.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-14-deep-iso0/img_avg_2015-11-18-14-deep-iso0_unit38_ch15_all_3sec_window_03893_spikes.npy'

#OK, but not great
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-7-9electrodein-iso0/img_avg_2015-11-18-7-9electrodein-iso0_unit07_ch16_modes_3sec_window_01191_spikes_first.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-7-9electrodein-iso0/img_avg_2015-11-18-7-9electrodein-iso0_unit07_ch16_modes_3sec_window_01191_spikes_tonic.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-7-9electrodein-iso0/img_avg_2015-11-18-7-9electrodein-iso0_unit07_ch16_modes_3sec_window_01191_spikes_last.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-7-9electrodein-iso0/img_avg_2015-11-18-7-9electrodein-iso0_unit07_ch16_modes_3sec_window_01191_spikes_burst.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-7-9electrodein-iso0/img_avg_2015-11-18-7-9electrodein-iso0_unit07_ch16_all_3sec_window_01191_spikes.npy'

#GREAT first vs. all
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-7-9electrodein-iso0/img_avg_2015-11-18-7-9electrodein-iso0_unit05_ch15_modes_3sec_window_04376_spikes_first.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-7-9electrodein-iso0/img_avg_2015-11-18-7-9electrodein-iso0_unit05_ch15_modes_3sec_window_04376_spikes_tonic.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-7-9electrodein-iso0/img_avg_2015-11-18-7-9electrodein-iso0_unit05_ch15_modes_3sec_window_04376_spikes_last.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-7-9electrodein-iso0/img_avg_2015-11-18-7-9electrodein-iso0_unit05_ch15_modes_3sec_window_04376_spikes_burst.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-7-9electrodein-iso0/img_avg_2015-11-18-7-9electrodein-iso0_unit05_ch15_all_3sec_window_04376_spikes.npy'

#GREAT first vs all; interesting bursting. 
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-7-9electrodein-iso0/img_avg_2015-11-18-7-9electrodein-iso0_unit01_ch11_modes_3sec_window_01345_spikes_first.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-7-9electrodein-iso0/img_avg_2015-11-18-7-9electrodein-iso0_unit01_ch11_modes_3sec_window_01345_spikes_tonic.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-7-9electrodein-iso0/img_avg_2015-11-18-7-9electrodein-iso0_unit01_ch11_modes_3sec_window_01345_spikes_last.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-7-9electrodein-iso0/img_avg_2015-11-18-7-9electrodein-iso0_unit01_ch11_modes_3sec_window_01345_spikes_burst.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-7-9electrodein-iso0/img_avg_2015-11-18-7-9electrodein-iso0_unit01_ch11_all_3sec_window_01345_spikes.npy'

#GREAT similar/identical modes
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-2/img_avg_2015-7-22-2_unit00_ch05_modes_3sec_window_00965_spikes_first.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-2/img_avg_2015-7-22-2_unit00_ch05_modes_3sec_window_00965_spikes_tonic.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-2/img_avg_2015-7-22-2_unit00_ch05_modes_3sec_window_00965_spikes_last.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-2/img_avg_2015-7-22-2_unit00_ch05_modes_3sec_window_00965_spikes_burst.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-2/img_avg_2015-7-22-2_unit00_ch05_all_3sec_window_00965_spikes.npy'

#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-2/img_avg_2015-7-22-2_unit04_ch15_modes_3sec_window_05350_spikes_first.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-2/img_avg_2015-7-22-2_unit04_ch15_modes_3sec_window_05350_spikes_tonic.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-2/img_avg_2015-7-22-2_unit04_ch15_modes_3sec_window_05350_spikes_last.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-2/img_avg_2015-7-22-2_unit04_ch15_modes_3sec_window_05350_spikes_burst.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-2/img_avg_2015-7-22-2_unit04_ch15_all_3sec_window_00794_spikes.npy'

#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-3/img_avg_2015-7-22-3_unit12_ch15_modes_3sec_window_06057_spikes_first.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-3/img_avg_2015-7-22-3_unit12_ch15_modes_3sec_window_06057_spikes_tonic.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-3/img_avg_2015-7-22-3_unit12_ch15_modes_3sec_window_06057_spikes_last.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-3/img_avg_2015-7-22-3_unit12_ch15_modes_3sec_window_06057_spikes_burst.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-3/img_avg_2015-7-22-3_unit12_ch15_all_3sec_window_06057_spikes.npy'


#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-1/img_avg_2015-7-22-1_unit01_ch06_modes_3sec_window_01690_spikes_first.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-1/img_avg_2015-7-22-1_unit01_ch06_modes_3sec_window_01690_spikes_tonic.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-1/img_avg_2015-7-22-1_unit01_ch06_modes_3sec_window_01690_spikes_last.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-1/img_avg_2015-7-22-1_unit01_ch06_modes_3sec_window_01690_spikes_burst.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-1/img_avg_2015-7-22-1_unit01_ch06_all_3sec_window_01690_spikes.npy'


#NOT GREAT
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-4/img_avg_2015-7-22-4_unit07_ch05_modes_3sec_window_00643_spikes_first.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-4/img_avg_2015-7-22-4_unit07_ch05_modes_3sec_window_00643_spikes_tonic.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-4/img_avg_2015-7-22-4_unit07_ch05_modes_3sec_window_00643_spikes_last.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-4/img_avg_2015-7-22-4_unit07_ch05_modes_3sec_window_00643_spikes_burst.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-4/img_avg_2015-7-22-4_unit07_ch05_all_3sec_window_00643_spikes.npy'



#GREAT MODE clustering; All modes identical!
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-5/img_avg_2015-7-22-5_unit15_ch12_modes_3sec_window_01152_spikes_first.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-5/img_avg_2015-7-22-5_unit15_ch12_modes_3sec_window_01152_spikes_tonic.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-5/img_avg_2015-7-22-5_unit15_ch12_modes_3sec_window_01152_spikes_last.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-5/img_avg_2015-7-22-5_unit15_ch12_modes_3sec_window_01152_spikes_burst.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-5/img_avg_2015-7-22-5_unit15_ch12_all_3sec_window_01152_spikes.npy'

#OK similar first and all
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-5/img_avg_2015-7-22-5_unit04_ch05_modes_3sec_window_01113_spikes_first.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-5/img_avg_2015-7-22-5_unit04_ch05_modes_3sec_window_01113_spikes_tonic.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-5/img_avg_2015-7-22-5_unit04_ch05_modes_3sec_window_01113_spikes_last.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-5/img_avg_2015-7-22-5_unit04_ch05_modes_3sec_window_01113_spikes_burst.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-5/img_avg_2015-7-22-5_unit04_ch05_all_3sec_window_01113_spikes.npy'

#NOt great
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-18-5electrodeinthalamus-iso0/img_avg_2015-12-1-18-5electrodeinthalamus-iso0_unit06_ch15_modes_3sec_window_01994_spikes_first.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-18-5electrodeinthalamus-iso0/img_avg_2015-12-1-18-5electrodeinthalamus-iso0_unit06_ch15_modes_3sec_window_01994_spikes_tonic.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-18-5electrodeinthalamus-iso0/img_avg_2015-12-1-18-5electrodeinthalamus-iso0_unit06_ch15_modes_3sec_window_01994_spikes_last.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-18-5electrodeinthalamus-iso0/img_avg_2015-12-1-18-5electrodeinthalamus-iso0_unit06_ch15_modes_3sec_window_01994_spikes_burst.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-18-5electrodeinthalamus-iso0/img_avg_2015-12-1-18-5electrodeinthalamus-iso0_unit06_ch15_all_3sec_window_01994_spikes.npy'

#OK, but modes are bit artificial
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-15-5electrodeinthalamus-iso0.8/img_avg_2015-12-1-15-5electrodeinthalamus-iso0.8_unit17_ch16_modes_3sec_window_02402_spikes_first.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-15-5electrodeinthalamus-iso0.8/img_avg_2015-12-1-15-5electrodeinthalamus-iso0.8_unit17_ch16_modes_3sec_window_02402_spikes_tonic.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-15-5electrodeinthalamus-iso0.8/img_avg_2015-12-1-15-5electrodeinthalamus-iso0.8_unit17_ch16_modes_3sec_window_02402_spikes_last.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-15-5electrodeinthalamus-iso0.8/img_avg_2015-12-1-15-5electrodeinthalamus-iso0.8_unit17_ch16_modes_3sec_window_02402_spikes_burst.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-15-5electrodeinthalamus-iso0.8/img_avg_2015-12-1-15-5electrodeinthalamus-iso0.8_unit17_ch16_all_3sec_window_02402_spikes.npy'


#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-15-5electrodeinthalamus-iso0.8/img_avg_2015-12-1-15-5electrodeinthalamus-iso0.8_unit18_ch16_modes_3sec_window_01527_spikes_first.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-15-5electrodeinthalamus-iso0.8/img_avg_2015-12-1-15-5electrodeinthalamus-iso0.8_unit18_ch16_modes_3sec_window_01527_spikes_tonic.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-15-5electrodeinthalamus-iso0.8/img_avg_2015-12-1-15-5electrodeinthalamus-iso0.8_unit18_ch16_modes_3sec_window_01527_spikes_last.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-15-5electrodeinthalamus-iso0.8/img_avg_2015-12-1-15-5electrodeinthalamus-iso0.8_unit18_ch16_modes_3sec_window_01527_spikes_burst.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-15-5electrodeinthalamus-iso0.8/img_avg_2015-12-1-15-5electrodeinthalamus-iso0.8_unit10_ch14_all_3sec_window_00279_spikes.npy'


#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-22-allelectrodeinthalamus-iso0/img_avg_2015-12-1-22-allelectrodeinthalamus-iso0_unit25_ch09_modes_3sec_window_04753_spikes_first.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-22-allelectrodeinthalamus-iso0/img_avg_2015-12-1-22-allelectrodeinthalamus-iso0_unit25_ch09_modes_3sec_window_04753_spikes_tonic.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-22-allelectrodeinthalamus-iso0/img_avg_2015-12-1-22-allelectrodeinthalamus-iso0_unit25_ch09_modes_3sec_window_04753_spikes_last.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-22-allelectrodeinthalamus-iso0/img_avg_2015-12-1-22-allelectrodeinthalamus-iso0_unit25_ch09_modes_3sec_window_04753_spikes_burst.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-22-allelectrodeinthalamus-iso0/img_avg_2015-12-1-22-allelectrodeinthalamus-iso0_unit25_ch09_all_3sec_window_00352_spikes.npy'


#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-22-allelectrodeinthalamus-iso0/img_avg_2015-12-1-22-allelectrodeinthalamus-iso0_unit27_ch10_modes_3sec_window_01742_spikes_first.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-22-allelectrodeinthalamus-iso0/img_avg_2015-12-1-22-allelectrodeinthalamus-iso0_unit27_ch10_modes_3sec_window_01742_spikes_tonic.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-22-allelectrodeinthalamus-iso0/img_avg_2015-12-1-22-allelectrodeinthalamus-iso0_unit27_ch10_modes_3sec_window_01742_spikes_last.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-22-allelectrodeinthalamus-iso0/img_avg_2015-12-1-22-allelectrodeinthalamus-iso0_unit27_ch10_modes_3sec_window_01742_spikes_burst.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-22-allelectrodeinthalamus-iso0/img_avg_2015-12-1-22-allelectrodeinthalamus-iso0_unit27_ch10_all_3sec_window_01742_spikes.npy'


#OK; used in data
'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-14-allelectrodeinthalamus-is0/img_avg_2015-12-2-14-allelectrodeinthalamus-is0_unit16_ch07_modes_3sec_window_00871_spikes_first.npy',
'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-14-allelectrodeinthalamus-is0/img_avg_2015-12-2-14-allelectrodeinthalamus-is0_unit16_ch07_modes_3sec_window_00871_spikes_tonic.npy',
'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-14-allelectrodeinthalamus-is0/img_avg_2015-12-2-14-allelectrodeinthalamus-is0_unit16_ch07_modes_3sec_window_00871_spikes_last.npy',
'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-14-allelectrodeinthalamus-is0/img_avg_2015-12-2-14-allelectrodeinthalamus-is0_unit16_ch07_modes_3sec_window_00871_spikes_burst.npy',
'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-14-allelectrodeinthalamus-is0/img_avg_2015-12-2-14-allelectrodeinthalamus-is0_unit16_ch07_all_3sec_window_00871_spikes.npy'


]



mid_line_mask = 10
block = 10



path_dir, fname = os.path.split(files[0])

unit = files[0][files[0].find("unit")+4:files[0].find("unit")+6]
print "... unit: ", unit

fname = glob.glob(path_dir+'/unit_'+str(unit).zfill(2)+"*_imagingspikes_grouped_spikes.txt")[0]

mode_spikes = []
with open(fname, 'rt') as inputfile:
    reader = csv.reader(inputfile)
    for row in reader:
        mode_spikes.append(np.float32(row))

mode_spikes.append(np.hstack(mode_spikes[:4]))

#****************************************************************8

clr_order = [0, 1, 2, 3, 4]

#spiking_modes = ['last', 'first', 'burst', 'tonic', 'all']
spiking_modes = ['first', 'last', 'burst', 'tonic', 'all']

colours = ['blue', 'green', 'red', 'black', 'magenta', 'orange','cyan','brown']
if '7-22-5' in files[0]: clr_order = [0, 1, 2, 3, 4]
if '7-22-2' in files[0]: clr_order = [3, 1, 0, 2, 4]
if '7-22-1' in files[0]: clr_order = [2, 1, 3, 0, 4]
if '12-2-14' in files[0]: clr_order = [1,0,2,3,4]

colours_temp = colours[clr_order[0]], colours[clr_order[1]], colours[clr_order[2]], colours[clr_order[3]], colours[4]
spiking_modes_temp = spiking_modes[clr_order[0]], spiking_modes[clr_order[1]], spiking_modes[clr_order[2]], spiking_modes[clr_order[3]], spiking_modes[4]
colours = colours_temp
spiking_modes = spiking_modes_temp
print colours
print spiking_modes


if False: 
    for k in range(len(files)):
        ax = plt.subplot(len(files),1,k+1)
        
        print "...group: ", k

        data = np.load(files[k])
        print data.shape

        temp_stack = []
        for p in range(0, len(data), block):
            temp_stack.append(np.mean(data[p:p+block], axis=0))

        temp_stack = quick_mask_single_allframe(temp_stack, mid_line_mask)
        
        img_out = np.ma.hstack((temp_stack))
        img_out[:, img_out.shape[1]/2-2:img_out.shape[1]/2+2] = np.min(img_out)
        
        v_max = np.nanmax(np.abs(img_out)); v_min = -v_max

        plt.imshow(img_out, vmin=v_min, vmax=v_max)
        
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        
        ax.set_title(spiking_modes[k] + " #spikes: " + str(len(mode_spikes[k])) + "  DF/F (max/min): " + str(round(v_max*100,1)) +"%", fontsize =30)

    plt.suptitle(files[0])
    plt.show()


#************************************** PLOT STM Time Dynamics ****************************************

area_names = ['barrel'] #, 'retrosplenial', 'motor', 'visual'] 
sides = ['left']#,'right']
depth = 'cortex'#, 'subcortical']

crop_method = 'max'  #,['max','ave']

dir_path, file_name = os.path.split(files[0]) 
vector_stack = []
ax=plt.subplot(111)
for ctr, file_ in enumerate(files):    
    side, area = sides[0], area_names[0]
    
    save_file = file_[:-4] + '_' + depth + "_" + area + "_" + side+"_roi"
    out_file = save_file+"_"+crop_method

    #Load data
    img_stack = np.load(file_)
    #ax=plt.subplot(131)
    #v_max = np.nanmax(np.abs(img_stack[0])); v_min = -v_max
    #plt.imshow(img_stack[0], vmin=-v_max, vmax=v_max)
    #print v_max
    
    #print img_stack.shape
    #temp_stack = []
    #for k in range(len(img_stack)):
    #    img_stack[k] = np.nan_to_num(img_stack[k])
    #    temp_stack.append(scipy.misc.imresize(img_stack[k], 64./len(img_stack[0])))      #Subsample data to fit 64 x 64 roi masks

    #ax=plt.subplot(132)
    #v_max = np.nanmax(np.abs(img_stack[0])); v_min = -v_max
    #plt.imshow(img_stack[0], vmin=-v_max, vmax=v_max)
    #print v_max

    #img_stack = np.array(temp_stack)
    #v_max = np.nanmax(np.abs(img_stack[0])); v_min = -v_max
    #ax=plt.subplot(133)
    #plt.imshow(img_stack[0], vmin=-v_max, vmax=v_max)
    #print v_max
    #plt.show()
        
    #Load mask
    mask_file = dir_path + '/' + depth + "_" + area + "_" + side       #Load ROI Mask
    temp_data = np.load(mask_file+'.npy')
    print temp_data.shape
    
    temp_data = scipy.misc.imresize(temp_data, (256, 256))      #Subsample data to fit 64 x 64 roi masks
    print temp_data.shape
    
    mask_stack_ravel = temp_data.ravel()
    print mask_stack_ravel.shape            #NB: MASK uses 64 x 64 subsampling; MUST CONVERT ALL DATA TO THIS
    
    indexes = np.where(mask_stack_ravel==0)[0]      #Pick only values not masked, i.e. mask = False
    
    #Loop over each frame and compute ave or max values in ROI
    area_vector_stack= []
    print "... area: ", area, " side: ", side
    area_vector = []
    for k in range(len(img_stack)):              #Loop over all frames 
        img_frame_ravel = img_stack[k].ravel()
        
        if crop_method == 'ave': area_vector.append(np.mean(img_frame_ravel[indexes]))
        if crop_method == 'max': area_vector.append(np.max(img_frame_ravel[indexes]))

    area_vector = np.array(area_vector)*100
    #Normalize the spike stack
    if False: 
        val_max = np.nanmax(area_vector); val_min = np.nanmin(area_vector)
        area_vector = (area_vector - val_min)/(val_max-val_min)


    plt.plot(area_vector, color = colours[ctr], linewidth = 5, alpha=0.8)


plt.plot([0,180],[0,0], 'r--', color='black', linewidth=3, alpha=0.5)
plt.plot([89,89],[-5,5], 'r--', color='black', linewidth=3, alpha=0.5)
plt.ylim(-1,1)
plt.xlim(0,178)

plt.ylabel("DF/F (%)", fontsize=30)

#plt.title("STMTD for Left-"+area+" Cortex", fontsize = 30)
plt.xlabel("Time (sec)", fontsize=30)

old_xlabel = np.arange(0, len(area_vector)+3,30)
new_xlabel = np.arange(-3.0,3.1, 1)
plt.xticks(old_xlabel, new_xlabel, fontsize=30)    

plt.tick_params(axis='both', which='both', labelsize=30)


#Plot legend
import matplotlib.patches as mpatches

blue = mpatches.Patch(facecolor = 'blue', edgecolor="black")
green = mpatches.Patch(facecolor = 'green', edgecolor="black")
red = mpatches.Patch(facecolor = 'red', edgecolor="black")
black = mpatches.Patch(facecolor = 'black', edgecolor="black")
magenta = mpatches.Patch(facecolor = 'magenta', edgecolor="black")

labels = spiking_modes

ax.legend([blue,green,red,black,magenta], labels, fontsize=12, loc=0)



plt.show()













