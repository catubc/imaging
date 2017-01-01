import numpy as np
import matplotlib.pyplot as plt
import struct
import string
import matplotlib.animation as animation

import matplotlib.gridspec as gridspec

#work_dir = '/media/cat/12TB/in_vivo/tim/allen/20150303/Evoked/L-HL_stim_1ms_1mA/'
#work_dir = '/media/cat/12TB/in_vivo/tim/allen/20150303/Evoked/L-FL stim_1ms 1mA/'
#work_dir = '/media/cat/12TB/in_vivo/tim/allen/20150303/Evoked/L-FL stim_1ms 2mA/'
#work_dir = '/media/cat/12TB/in_vivo/tim/allen/20150303/Evoked/L-whisker stim_1ms 5V/'
#work_dir = '/media/cat/12TB/in_vivo/tim/allen/20150303/Evoked/L-visual stim_1ms/'
#work_dir = '/media/cat/12TB/in_vivo/tim/allen/20150303/Evoked/L-aud stim_1ms/'
#file_names = ['hi1', 'hi2', 'hi3', 'hi4', 'hi5','lo1', 'lo2', 'lo3', 'lo4', 'lo5']

file_names = [
'/media/cat/12TB/in_vivo/tim/yuki/AQ2/AQ2am_Jan18_30Hz_reward_code_07.npy',
'/media/cat/12TB/in_vivo/tim/yuki/AQ2/AQ2am_Jan18_30Hz_reward_code_04.npy',
'/media/cat/12TB/in_vivo/tim/yuki/AQ2/AQ2am_Jan18_30Hz_reward_code_02.npy'
]

#Load General mask (removes background)
if False:
    n_pixels = 128
    generic_mask_file = []
    generic_mask_file = work_dir + 'genericmask.txt'
    if (os.path.exists(generic_mask_file)==True):
        generic_coords = np.loadtxt(generic_mask_file)
    else:
        generic_coords = Define_generic_mask(data, work_dir, file_name, img_rate, n_pixels)

    generic_mask_indexes=np.zeros((n_pixels,n_pixels))
    for i in range(len(generic_coords)):
        generic_mask_indexes[generic_coords[i][0]][generic_coords[i][1]] = True

    data = []
    for i in range(len(file_names)):
        data_load = np.load(work_dir+file_names[i]+'_processed.npy')
        v_min = np.min(data_load)
        v_max = np.max(data_load)
        data_load = (data_load - v_min)/(v_max-v_min)

        temp_array = []
        #Mask all frames; NB: PROBABLY FASTER METHOD
        for i in range(0, len(data_load),1):
            temp = np.ma.array(data_load[i], mask=generic_mask_indexes, fill_value = 0, hard_mask = True)
            temp_array.append(temp)

        data.append(temp_array)

#Load data
data = []
for k in range(len(file_names)):
    data.append(np.load(file_names[k]))



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
        
