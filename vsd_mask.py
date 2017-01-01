import matplotlib.pyplot as plt
import numpy as np
import os


files = [

'/media/cat/12TB/in_vivo/tim/cat/2016_07_11_vsd/tif_files/track1_150Hz_iso1.0_spontaneous.npy',
#'/media/cat/12TB/in_vivo/tim/cat/2016_07_11_vsd/tif_files/track2_150Hz_iso1.0_spontaneous.npy',

'/media/cat/12TB/in_vivo/tim/cat/2016_07_12_vsd/tif_files/track1_150Hz_1_160712_202358.npy',
#'/media/cat/12TB/in_vivo/tim/cat/2016_07_12_vsd/tif_files/track2_150Hz_1_160712_213342.npy',

'/media/cat/12TB/in_vivo/tim/cat/2016_07_15_vsd/tif_files/track1_fsf_200_repeats_160715_204633.npy',
#'/media/cat/12TB/in_vivo/tim/cat/2016_07_15_vsd/tif_files/track1_150Hz_1st_spontaneous_10iso_160715_181445.npy',
#'/media/cat/12TB/in_vivo/tim/cat/2016_07_15_vsd/tif_files/track1_150Hz_2nd_spontaneous_15iso_160715_192621.npy',
#'/media/cat/12TB/in_vivo/tim/cat/2016_07_15_vsd/tif_files/track1_150Hz_3rd_spontaneous_160715_211403.npy',
#'/media/cat/12TB/in_vivo/tim/cat/2016_07_15_vsd/tif_files/track2_spontaneous_1_160715_222854.npy',


#'/media/cat/500GB/in_vivo/tim/cat/2016_07_20_vsd/tif_files/track1_iso15_spontaneous_160720_210136.npy',
#'/media/cat/500GB/in_vivo/tim/cat/2016_07_20_vsd/tif_files/track2_iso75_spontaneous_2nd_160720_231526.npy',
#'/media/cat/500GB/in_vivo/tim/cat/2016_07_20_vsd/tif_files/track2_iso75_spontaneous_160720_220333.npy'

]

for file_ in files:

    if (os.path.exists(file_[:-4]+'_mask.npy')==False):
        print "...reading file: ", file_
        images_raw=0.; images_raw = np.load(file_) #.astype(np.float32, copy=False)
        #images_rawframe = np.average(images_raw[np.random.randint(0,len(images_raw),20000)], axis=0)    #Average 20000 frames for baseline
        print images_raw.shape
        if len(images_raw.shape)<3:
            print "..reshaping..."
            images_raw = images_raw.reshape(-1, 128, 128)
            print images_raw.shape
            
        print "...computing std..."
        images_std = np.std(images_raw, axis=0)
        plt.imshow(images_std, vmin=np.min(images_std),vmax=np.max(images_std))
        plt.show()
    
        print "...saving to file..."
        np.save(file_[:-4]+'_std', images_std)
    
