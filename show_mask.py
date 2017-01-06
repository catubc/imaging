from matplotlib import pyplot as plt
import numpy as np
import scipy
import scipy.misc


data = np.zeros((64,64), dtype=np.float32)
midline_mask_n_pixels = 2
n_pixels = len(data[0])
    
generic_mask_file = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-11/genericmask.txt'
generic_coords = np.int32(np.loadtxt(generic_mask_file))

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
temp_array = np.ma.array(np.zeros((n_pixels,n_pixels),dtype=np.float32), mask=True)
temp_array = np.ma.masked_array(data, mask=generic_mask_indexes, fill_value=0)

plt.imshow(temp_array)
plt.show()

base_filename= '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-15-allelectrodeinthalamus-is0/cortex_'

areas=['motor', 'barrel', 'retrosplenial']
side = 'left'

for area in areas:
    mask = np.load(base_filename + area+"_"+side+".npy")

    print mask.shape
    for k in range(len(mask)):
        indexes = np.where(mask[k]==0)[0]
        temp_array[k][indexes]=1
        
plt.imshow(temp_array, interpolation='none')
plt.show()
