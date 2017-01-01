import matplotlib.pyplot as plt
import numpy as np
import os

files = [

#'/media/cat/12TB/in_vivo/tim/cat/2016_07_11_vsd/tif_files/track1_150Hz_iso1.0_spontaneous.npy',
#'/media/cat/12TB/in_vivo/tim/cat/2016_07_11_vsd/tif_files/track2_150Hz_iso1.0_spontaneous.npy',

#'/media/cat/12TB/in_vivo/tim/cat/2016_07_12_vsd/tif_files/track1_150Hz_1_160712_202358.npy',
#'/media/cat/12TB/in_vivo/tim/cat/2016_07_12_vsd/tif_files/track2_150Hz_1_160712_213342.npy',

##'/media/cat/12TB/in_vivo/tim/cat/2016_07_15_vsd/tif_files/track1_fsf_200_repeats_160715_204633.npy',
#'/media/cat/12TB/in_vivo/tim/cat/2016_07_15_vsd/tif_files/track1_150Hz_1st_spontaneous_10iso_160715_181445.npy',
#'/media/cat/12TB/in_vivo/tim/cat/2016_07_15_vsd/tif_files/track1_150Hz_2nd_spontaneous_15iso_160715_192621.npy',
#'/media/cat/12TB/in_vivo/tim/cat/2016_07_15_vsd/tif_files/track1_150Hz_3rd_spontaneous_160715_211403.npy',
#'/media/cat/12TB/in_vivo/tim/cat/2016_07_15_vsd/tif_files/track2_spontaneous_1_160715_222854.npy',

#'/media/cat/500GB/in_vivo/tim/cat/2016_07_20_vsd/tif_files/track1_iso15_spontaneous_160720_210136.npy',
'/media/cat/500GB/in_vivo/tim/cat/2016_07_20_vsd/tif_files/track2_iso75_spontaneous_2nd_160720_231526.npy',
'/media/cat/500GB/in_vivo/tim/cat/2016_07_20_vsd/tif_files/track2_iso75_spontaneous_160720_220333.npy'

]

for file_ in files:

    if True:
    #if (os.path.exists(file_[:-4]+'_filtered.npy')==False):
        print "...reading file: ", file_
        images_raw=0.; images_raw = np.load(file_).astype(np.float64, copy=False)
        #images_rawframe = np.average(images_raw[np.random.randint(0,len(images_raw),20000)], axis=0)    #Average 20000 frames for baseline
        print "... making average..."
        images_rawframe = np.average(images_raw, axis=0)

        print "Temporal filtering ..."
        from scipy.signal import butter, filtfilt

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

        lowcut = 0.5
        highcut=6
        img_rate = 150.64        #Frame rate of imaging

        print "...band pass filtering: ", lowcut, "hz to ", highcut, "hz"
        images_filtered = np.zeros(images_raw.shape, dtype=np.float32)
        for i in range(len(images_filtered[0])):
            for j in range(len(images_filtered[0][0])):
                print "...filtering: ", i,j
                images_filtered[:, i, j] = np.array(butter_bandpass_filter(images_raw[:,i,j], lowcut, highcut, img_rate, order = 2))
        
        #n_splits = 20
        #for k in range(n_splits):
        #    print "...filtering chunk: ", k
        #    images_raw[k*len(images_raw)/n_splits:(k+1)*len(images_raw)/n_splits] = \
        #    butter_bandpass_filter(images_raw[k*len(images_raw)/n_splits:(k+1)*len(images_raw)/n_splits], lowcut, highcut, img_rate, order = 2)
      
        print "... adding base frame..."
        images_filtered += images_rawframe
    
        print "...saving to file..."
        np.save(file_[:-4]+'_filtered', images_filtered.astype(np.float32, copy=False))
    
