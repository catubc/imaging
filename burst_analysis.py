import matplotlib.pyplot as plt
import numpy as np
import os, csv
from sklearn import cluster, datasets

colours = ['blue', 'green', 'red', 'black', 'magenta', 'orange','cyan','brown']
for k in range(100): colours.append('brown')
[
    [230, 159,   0], #  41 orange
    [213,  94,   0], #  27 vermillion
    [  0, 114, 178], # 202 blue
    [204, 121, 167], # 326 reddish purple
    [  0,   0,   0], # --- black
    [ 86, 180, 233], # 202 sky blue
    [  0, 158, 115], # 164 bluish green
    [240, 228,  66], #  56 yellow
]

colours = np.array(colours)

#spiking_modes = {'burst':'blue', 'last':'green', 'tonic':'red', 'first':'black']
spiking_modes = ['first', 'last', 'burst', 'tonic']
print colours[:4]
print spiking_modes
#
quantile = 0.15

clr_order = [0,1,2,3] #Default order unless otherwise indicated

#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-8-9electrodein-iso0/unit_02_channel_11_ptp_149_imagingspikes.txt'          ; clr_order = [0,1,2,3]

#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-2/unit_00_channel_05_ptp_055_imagingspikes.txt'                               ; clr_order = [3,1,0,2]
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-2/unit_04_channel_15_ptp_300_imagingspikes.txt'                               ; clr_order = [0,3,2,1]
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-2/unit_02_channel_06_ptp_042_imagingspikes.txt'
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-2/unit_01_channel_05_ptp_062_imagingspikes.txt'

#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-3/unit_13_channel_15_ptp_159_imagingspikes.txt'
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-3/unit_12_channel_15_ptp_268_imagingspikes.txt'
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-3/unit_16_channel_16_ptp_093_imagingspikes.txt'

#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-1/unit_01_channel_06_ptp_037_imagingspikes.txt'                   ; clr_order = [2,1,3,0]

#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-5-10electrodein-iso0/unit_12_channel_16_ptp_053_imagingspikes.txt'         ; clr_order = [0,1,2,3]

#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-15-allelectrodeinthalamus-is0/unit_33_channel_11_ptp_078_imagingspikes.txt'  ; clr_order = [0,1,2,3]
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-12-allelectrodeinthalamus-iso1/unit_65_channel_16_ptp_074_imagingspikes.txt' ; clr_order = [0,1,2,3]
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-11/2015-12-11-12-allelectinthalamus-iso1/unit_02_channel_03_ptp_088_imagingspikes.txt'   ; clr_order = [0,1,2,3]
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-14-deep-iso0/unit_16_channel_09_ptp_140.csv'                                

#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-8-9electrodein-iso0/unit_03_channel_13_ptp_199_imagingspikes.txt'

#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/unit_01_channel_07_ptp_074_imagingspikes.txt'
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/unit_02_channel_08_ptp_084_imagingspikes.txt'
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/unit_00_channel_06_ptp_065_imagingspikes.txt'
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/unit_08_channel_11_ptp_077_imagingspikes.txt'
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-16-deep-iso0/unit_14_channel_12_ptp_065_imagingspikes.txt'

#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-13-deep-iso0/unit_13_channel_09_ptp_073_ds__imagingspikes.txt'
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-13-deep-iso0/unit_14_channel_10_ptp_064_imagingspikes.txt'

#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-14-deep-iso0/unit_07_channel_06_ptp_157_imagingspikes.txt'
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-14-deep-iso0/unit_09_channel_07_ptp_063_imagingspikes.txt'
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-14-deep-iso0/unit_21_channel_10_ptp_086_imagingspikes.txt'
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-14-deep-iso0/unit_35_channel_14_ptp_257_imagingspikes.txt'
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-14-deep-iso0/unit_37_channel_15_ptp_383_imagingspikes.txt'
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-14-deep-iso0/unit_38_channel_15_ptp_051_imagingspikes.txt'

#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-7-9electrodein-iso0/unit_07_channel_16_ptp_082_imagingspikes.txt'
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-7-9electrodein-iso0/unit_05_channel_15_ptp_171_imagingspikes.txt'
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-7-9electrodein-iso0/unit_01_channel_11_ptp_183_imagingspikes.txt'

#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-4/unit_07_channel_05_ptp_101_imagingspikes.txt'

#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-5/unit_15_channel_12_ptp_060_imagingspikes.txt'       ; clr_order = [0, 1,2,3]
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-5/unit_04_channel_05_ptp_093_imagingspikes.txt'

#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-18-5electrodeinthalamus-iso0/unit_06_channel_15_ptp_074_imagingspikes.txt'    ; clr_order = [0,1,2,3]
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-18-5electrodeinthalamus-iso0/unit_09_channel_16_ptp_068_imagingspikes.txt'

#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-15-5electrodeinthalamus-iso0.8/unit_17_channel_16_ptp_057_imagingspikes.txt'
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-15-5electrodeinthalamus-iso0.8/unit_18_channel_16_ptp_076_imagingspikes.txt'

#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-22-allelectrodeinthalamus-iso0/unit_17_channel_07_ptp_072_imagingspikes.txt'
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-22-allelectrodeinthalamus-iso0/unit_25_channel_09_ptp_126_imagingspikes.txt'
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-22-allelectrodeinthalamus-iso0/unit_26_channel_09_ptp_078_imagingspikes.txt'
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-1/2015-12-1-22-allelectrodeinthalamus-iso0/unit_27_channel_10_ptp_071_imagingspikes.txt'


#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-14-allelectrodeinthalamus-is0/unit_10_channel_05_ptp_074_imagingspikes.txt'
#file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-14-allelectrodeinthalamus-is0/unit_11_channel_06_ptp_075_imagingspikes.txt'
file_name_all = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-14-allelectrodeinthalamus-is0/unit_16_channel_07_ptp_238_imagingspikes.txt';  clr_order = [1,0,2,3]

colours_temp = colours[clr_order[0]], colours[clr_order[1]], colours[clr_order[2]], colours[clr_order[3]]
spiking_modes_temp = spiking_modes[clr_order[0]], spiking_modes[clr_order[1]], spiking_modes[clr_order[2]], spiking_modes[clr_order[3]]
colours = colours_temp
spiking_modes = spiking_modes_temp
print colours
print spiking_modes

spikes = np.loadtxt (file_name_all)

path_dir, filename = os.path.split(file_name_all) 


ax = plt.subplot(111)

temp_isi = []
temp_isi2 = []
f_rate=[]
raster_temp = spikes #Select raster from the state/depth separated and class separated lists
for q in range(1,len(raster_temp)-2,1):
    pre_spike = (raster_temp[q]-raster_temp[q-1])
    post_spike = (raster_temp[q+1]-raster_temp[q])
    if pre_spike == 0: continue     #Duplicates in data
    if post_spike == 0: continue    #Duplicates in data
    #if pre_spike < 0. or post_spike <0: print pre_spike, post_spike; quit()
    temp_isi.append(pre_spike)
    temp_isi2.append(post_spike)


cluster_data = np.log10(np.vstack((temp_isi,temp_isi2))).T
print cluster_data


val_thrshold = -0.5
if False: #Manual labeling
    labels = []
    for k in range(len(cluster_data)):
        if (cluster_data[k][0]<val_thrshold) and (cluster_data[k][1]<val_thrshold):
            labels.append(0)
        if (cluster_data[k][0]<val_thrshold) and (cluster_data[k][1]>val_thrshold):
            labels.append(1)
        if (cluster_data[k][0]>val_thrshold) and (cluster_data[k][1]>val_thrshold):
            labels.append(2)
        if (cluster_data[k][0]>val_thrshold) and (cluster_data[k][1]<val_thrshold):
            labels.append(3)
            
            
if True: #KMEANS
    n_clusters = 4
    clusters = cluster.KMeans(n_clusters, max_iter=1000, n_jobs=-1, random_state=1032)
    clusters.fit(cluster_data)

    labels = clusters.labels_

if False: #MEAN SHIFT
    from sklearn.cluster import MeanShift, estimate_bandwidth
    from sklearn.datasets.samples_generator import make_blobs
    
    quantile = 0.2
    bandwidth = estimate_bandwidth(cluster_data, quantile=quantile, n_samples=5000)

    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    ms.fit(cluster_data)
    labels = ms.labels_
    #print labels

if False: 
    from sklearn.cluster import DBSCAN
    from sklearn import metrics
    from sklearn.datasets.samples_generator import make_blobs
    from sklearn.preprocessing import StandardScaler 

    X = StandardScaler().fit_transform(cluster_data)

    eps = 0.2
    
    db = DBSCAN(eps=eps, min_samples=10).fit(X)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_ 

    #updated_data = []
    #for k in range(len(cluster_data)):
        #if labels[k]!=0:
            #updated_data.append(cluster_data[k])
    
    #cluster_data = np.array(updated_data)
    #X = StandardScaler().fit_transform(cluster_data)

    #db = DBSCAN(eps=0.1, min_samples=10).fit(X)
    #core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    #core_samples_mask[db.core_sample_indices_] = True
    #labels = db.labels_ 



labels = np.array(labels)

with open(file_name_all[:-4]+"_grouped_spikes.txt", "wt") as f:
    writer = csv.writer(f)
    for k in range(4): 
        print "... spiking mode: ", spiking_modes[k], "  #spikes: ", len(np.where(labels==k)[0])
        writer.writerows([spikes[np.where(labels==k)[0]]])
    

clrs = []
for k in range(len(labels)):
    clrs.append(colours[labels[k]])


plt.scatter(np.power(10,cluster_data[:,0]), np.power(10, cluster_data[:,1]), s=20, color =clrs)

plt.plot([1E-3, 1E1],[1E-1,1E-1], 'r--', color='black', linewidth =3, alpha=0.5)
plt.plot([1E-1,1E-1],[1E-3, 1E1], 'r--', color='black', linewidth =3, alpha=0.5)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1E-3,1E1)
ax.set_ylim(1E-3,1E1)
plt.suptitle(file_name_all)

plt.tick_params(axis='both', which='both', labelsize = 30)
plt.xlabel("ISI - Previous Spike (sec)", fontsize = 30)
plt.ylabel("ISI - Next Spike (sec)", fontsize = 30)
plt.show()
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
