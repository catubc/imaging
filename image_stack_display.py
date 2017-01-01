import numpy as np
import matplotlib.pyplot as plt
import struct
import string
import matplotlib.animation as animation
import os

import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable #Used for colorbar; allocates a bit of space for it

#work_dir = '/media/cat/12TB/in_vivo/tim/allen/20150303/Evoked/L-HL_stim_1ms_1mA/'
#file_names = ['hi1', 'hi2', 'hi3', 'hi4', 'hi5','lo1', 'lo2', 'lo3', 'lo4', 'lo5']

file_names = [
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-8-9electrodein-iso0/img_avg_2015-11-18-8-9electrodein-iso0_unit02_ch11_all_3sec_window_00453_spikes.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-2/img_avg_2015-7-22-2_unit00_ch05_all_3sec_window_00965_spikes.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-8-9electrodein-iso0/img_avg_2015-11-18-8-9electrodein-iso0_unit00_ch10_all_3sec_window_00704_spikes.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-5-10electrodein-iso0/img_avg_2015-11-27-5-10electrodein-iso0_unit12_ch16_all_3sec_window_01888_spikes.npy'
'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-5-10electrodein-iso0/img_avg_2015-11-27-5-10electrodein-iso0_unit12_ch16_all_3sec_window_01888_spikes.npy'
##'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-7-10electrodeincortex-iso0/img_avg_2015-12-2-7-10electrodeincortex-iso0_unit10_ch15_all_3sec_window_06101_spikes.npy',


###subcortical cells:
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-15-allelectrodeinthalamus-is0/img_avg_2015-12-2-15-allelectrodeinthalamus-is0_unit33_ch11_all_3sec_window_04290_spikes.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-12-allelectrodeinthalamus-iso1/img_avg_2015-12-2-12-allelectrodeinthalamus-iso1_unit65_ch16_all_3sec_window_00448_spikes.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-11/2015-12-11-12-allelectinthalamus-iso1/img_avg_2015-12-11-12-allelectinthalamus-iso1_unit02_ch03_all_3sec_window_01213_spikes.npy',
#'/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-14-deep-iso0/img_avg_2015-11-18-14-deep-iso0_unit16_ch09_all_3sec_window_02036_spikes.npy',

]

#Load General mask (removes background)
if True:
    n_pixels = 256
    generic_mask_file = []
    generic_mask_file = '/media/cat/8TB/in_vivo/tim/dongsheng/genericmask.txt'
    generic_coords = np.loadtxt(generic_mask_file)
    print len(generic_coords)

    generic_mask_indexes=np.zeros((n_pixels,n_pixels))
    for i in range(len(generic_coords)):
        generic_mask_indexes[generic_coords[i][0]][generic_coords[i][1]] = True
            
data = []
for ctr, file_name in enumerate(file_names):
    ax = plt.subplot(4,1, ctr+1)
    
    data_load = np.load(file_name)

    temp_array = []
    #Mask all frames; NB: PROBABLY FASTER METHOD
    for i in range(0, len(data_load),1):
        temp = np.ma.array(data_load[i], mask=generic_mask_indexes, fill_value = 0, hard_mask = True)
        temp_array.append(temp)
    data = temp_array

    #v_max= np.nanmax(data);  v_min = np.nanmin(data); 
    #if v_min>-v_max: v_max = -v_min
    #else: v_min = -v_max
    #print v_min, v_max

    block = 17
    temp_stack = []
    for q in range(0,len(data), block):
        temp_stack.append(np.ma.average(data[q:q+block], axis=0))
    temp_plot = np.ma.hstack(temp_stack); v_max= np.max(np.abs(temp_plot));  v_min = -v_max
    print v_min, v_max

    ax.yaxis.set_ticks([])
    right_border = 3.0; left_border = -3.0
    old_xlabel = np.linspace(0,temp_plot.shape[1], (right_border - left_border)*10/block+1)
    new_xlabel = np.around(np.linspace(left_border,right_border, (right_border - left_border)*10/block+1), decimals=2)
    plt.xticks(old_xlabel, new_xlabel, fontsize=15)

    temp_plot[:][len(temp_plot[0])/2-3:len(temp_plot[0])/2+3]=v_min
    im = plt.imshow(temp_plot, vmin=v_min, vmax=v_max)
    plt.title(file_name[file_name.index('img'):])

    #Colorbar
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", "2%", pad="1%")
    cbar = plt.colorbar(im, cax=cax, ticks=[v_min, 0, v_max], orientation='vertical')
    cbar.ax.set_yticklabels([str(round(v_min*100,1))+'%', '0%', str(round(v_max*100,1))+'%'], fontsize=15)  # horizontal colorbar
    
    
plt.show()


#Data stack;
gs = gridspec.GridSpec(2, 2)
ax1 = plt.subplot(gs[0,0])
plt.imshow(data[0][0])

ax1 = plt.subplot(gs[0,1])
plt.imshow(data[1][0])

ax1 = plt.subplot(gs[1,:])
test1 = np.ma.dstack((data[0],data[1]))
plt.imshow(test1[0])
plt.show()

data = np.dstack((data[0],data[1],data[2]))
#row2 = np.ma.dstack((data[5],data[6],data[7],data[8],data[9]))

#data = np.ma.hstack((row1,row2))

v_max=np.max(data)
v_min=np.min(data)

Writer = animation.writers['ffmpeg']
writer = Writer(fps=5, metadata=dict(artist='Me'), bitrate=1800)

fig = plt.figure() # make figure

# make axesimage object
# the vmin and vmax here are very important to get the color map correct
im = plt.imshow(data[0], cmap=plt.get_cmap('jet'), interpolation='none')#, vmin=0, vmax=v_max)
# function to update figure
img_rate = 30.6
window = int(2*img_rate)
def updatefig(j):
    # set the data in the axesimage object
    im.set_array(data[j])
    #if 'visual' in file_names[0]:
    #    plt.title("Frame: "+str(j+34)+"\n"+str(round(float(j+34)/150,2))+"sec")
    #else:
    #    plt.title("Frame: "+str(j)+"\n"+str(round(float(j)/150,2))+"sec")

    plt.title("Frame: "+str(j)+"  " +str(round((float(j)-window)/img_rate,2))+"sec")

    # return the artists set
    return im,
# kick off the animation
ani = animation.FuncAnimation(fig, updatefig, frames=range(len(data)), interval=100, blit=False, repeat=True)

ani.save(file_names[0]+'_combined.mp4', writer=writer)

plt.show()
        
