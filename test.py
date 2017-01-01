import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import scipy.ndimage as ndimage

#*******************************
def make_vids(data, file_, main_dir):

    for k in range(len(data[0])):
        data[0][k]= ndimage.gaussian_filter(data[0][k], sigma=1)

    data = mask_data(main_dir, filename, data)

    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=6400)

    fig = plt.figure() # make figure

    im = []
    for k in range(len(data)):
        ax = plt.subplot(1,1,k+1)

        ax.get_xaxis().set_visible(False)
        ax.yaxis.set_ticks([])
        ax.yaxis.labelpad = 0
        Ev_max = np.max(np.abs(data[k]))*.5
        #v_min = -v_max 
        print "Max/mind DFF (%): ", Ev_max*100
        
        im.append([])
        im[k] = plt.imshow(data[k][0], cmap=plt.get_cmap('jet'), vmin = -Ev_max, vmax = Ev_max, interpolation='none')#, vmin=0, vmax=v_max)

    #function to update figure
    def updatefig(j):
        if (j%100==0): print "...frame: ", j
        plt.suptitle(filename.replace('/media/cat/500GB/in_vivo/tim/cat/2016_07_20_vsd/stm_files/','')+
        "\nFrame: "+str(j)+"\n"+str(round(float(j)/len(data[0])*(len(data[0])/150.64)-len(data[0])/2./150.64,2))+"sec", fontsize=10)

        # set the data in the axesimage object
        for k in range(len(data)):
            im[k].set_array(data[k][j])

        # return the artists set
        return im
        
    ani = animation.FuncAnimation(fig, updatefig, frames=range(len(data[0])), interval=10, blit=False, repeat=True)
    ani.save(file_+'.mp4', writer=writer)
    plt.show()

#********************************
def mask_data(main_dir, file_, data):
    print "Masking data..."
    img_rate = 150      #Frame rate of imaging
    n_pixels = len(data[0][0])      #number of pixels

    generic_mask_file = []
    generic_mask_file = main_dir + 'genericmask.txt'
    generic_coords = np.loadtxt(generic_mask_file)

    generic_mask_indexes=np.zeros((n_pixels,n_pixels))
    for i in range(len(generic_coords)):
        generic_mask_indexes[int(generic_coords[i][0])][int(generic_coords[i][1])] = True

    for k in range(len(data)):
        temp_array = []
        for i in range(0, len(data[k]),1):
            temp = np.ma.array(data[k][i], mask=generic_mask_indexes, fill_value = 0., hard_mask = True)
            temp_array.append(temp)
        data[k] = np.ma.array(temp_array)#[:100]
    
    return data
    
#***********************
main_dir = '/media/cat/12TB/in_vivo/tim/cat/2016_07_11_vsd/'

#Load General mask (removes background)
generic_mask_file = []
generic_mask_file = main_dir + 'genericmask.txt'
generic_coords = np.int16(np.loadtxt(generic_mask_file))

generic_mask_indexes=np.zeros((128,128))
for i in range(len(generic_coords)):
    generic_mask_indexes[generic_coords[i][0]][generic_coords[i][1]] = True

#*******************************
filename = '/media/cat/12TB/in_vivo/tim/cat/2016_07_11_vsd/stm_files/img_avg_track1_150Hz_iso1.0_spontaneous_lp_compressed_unit001_ch028_all_3sec_window_00147_spikes.npy'

print "Loading file: ", filename
images_filtered = np.load(filename)
print images_filtered.shape


#print np.isnan(images_filtered).any()
#quit()


print "...computing baseline..."
data=[]
#images_filtered1 = images_filtered - images_filtered[3000]
#baseline = np.average(images_filtered[300:], axis=0)       #Cat's method: uses minimum offset for processing
#data1 = (images_filtered[300:]-baseline)/baseline
#data.append(data1[:][:1000])
data.append(images_filtered)
#data = mask_data(main_dir, filename, data)


make_vids([data[0][:]], filename, main_dir)
quit()

plt.imshow(np.ma.array(img1[450],mask=generic_mask_indexes))
plt.show()

for k in range(0,9,1):
    ax = plt.subplot(3,3,k+1)
    plt.imshow(np.ma.array(img1[k*10+400],mask=generic_mask_indexes))
    #plt.imshow(img1[400+k])
plt.show()
quit()

img2 = img1



for k in range(len(img1)):
    frame = k
    ax=plt.subplot(121)
    v_min = np.nanmin(img1[frame])
    v_max = np.nanmax(img1[frame])
    ax.imshow(np.ma.array(img1[frame],mask=generic_mask_indexes), vmin=v_min, vmax=v_max)

    ax=plt.subplot(122)
    v_min = np.nanmin(img2[frame])
    v_max = np.nanmax(img2[frame])
    ax.imshow(np.ma.array(img2[frame],mask=generic_mask_indexes), vmin=v_min, vmax=v_max)
    plt.show()



