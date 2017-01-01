import numpy as np
import matplotlib.pyplot as plt
import glob


base_filename = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-2/stack1D_2015-7-22-2_unit'



cells = np.arange(0,100,1)
cells = [2]

methods_array = np.arange(0,3,1)
methods_array = [0]

#partition_data(base_filename, cells, methods, methods_array, colors)

filename = glob.glob(base_filename+str(cells[0]).zfill(2)+"*spikes.npy")[0]


original_image1 = np.load(filename, mmap_mode='r')

plt.imshow(original_image1[0,:,:64])

plt.show()
