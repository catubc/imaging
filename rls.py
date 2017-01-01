import numpy as np
import matplotlib.pylab as plt
import padasip as pa 
import os
import scipy, scipy.misc

#************************************************************************************

def quick_mask_single_frame(data, midline_mask_n_pixels):
    
    n_pixels = len(data[0])
        
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
    temp_array = np.ma.array(np.zeros((n_pixels,n_pixels),dtype=np.float32), mask=True)
    temp_array = np.ma.masked_array(data, mask=generic_mask_indexes, fill_value=0)
    
    return temp_array


# these two function supplement your online measurment
def measure_x():
    # it produces input vector of size 3
    x = np.random.random(3)
    return x

def measure_d(x, y_total):
    # meausure system output
    #d = 2*x[0] + 1*x[1] - 1.5*x[2]
    
    d = x - y_total

    return d
    
    
    
#************************************************************************************

filename = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-2/stack1D_2015-7-22-2_unit01_ch05_all_3sec_window_01141_spikes.npy'
#filename = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-6-10electrodeincortex-iso0/stack1D_2015-12-2-6-10electrodeincortex-iso0_unit09_ch14_all_3sec_window_01728_spikes.npy'

raw_data = np.load(filename, mmap_mode='c')
print "... # frames: ", raw_data.shape[2]/64
for f in range(raw_data.shape[2]/64):
    original_image = raw_data[:, :, f*64:(f+1)*64]

    print "...processing frame: ", f

    N = original_image.shape[0]
    log_d = np.zeros(N)
    log_y = np.zeros(N)
    y_total = np.zeros(len(original_image[0].ravel()))
    filter_array = []

    filt = pa.filters.FilterRLS(len(original_image[0].ravel()), mu=0.5)

    save_file = filename[:-4]+"_frame_"+str(f)+"_rls_array"
    
    if os.path.exists(save_file+'.npy')==False: 

        #ax = plt.subplots(1,1)
        for k in range(N):
            print "...spike: ", k
            # measure input
            #x = measure_x()
            x = original_image[k].ravel()
            
            # predict new value
            #print "...predicting new val: "
            y = filt.predict(x)
            #print y
            # do the important stuff with prediction output
            pass    
            
            # measure output
            #keep track of running total
            y_total = y_total + x*y

            d = measure_d(x, y_total)
            #d = y - x
            
            # update filter
            #print "...adapting filter: "
            filt.adapt(d, x)
            
            #log values
            #log_d[k] = d
            #log_y[k] = y
            
            #plt.clf()

            #plt.subplot(1,2,1)
            #plt.imshow(original_image[k])
            #ax=plt.subplot(1,2,2)

            temp_y = np.vstack(np.split(y_total,64))
            filter_array.append(temp_y)
            
            #temp_y = quick_mask_single_frame(temp_y, 1)

            #image.set_data(temp_y)

            # pp.draw()
            #fig.canvas.draw()

            #if k%5==0:
                #plt.imshow(temp_y)
                #plt.show()

        np.save(save_file, filter_array)
        #np.save(filename[:-4]+"_frame_"+str(f)+"_log_d_array", log_d)
        #np.save(filename[:-4]+"_frame_"+str(f)+"_log_y_array", log_y)



        
### show results
plt.figure(figsize=(15,9))
plt.subplot(211);plt.title("Adaptation");plt.xlabel("samples - k")
plt.plot(log_d,"b", label="d - target")
plt.plot(log_y,"g", label="y - output");plt.legend()
plt.subplot(212);plt.title("Filter error");plt.xlabel("samples - k")
plt.plot(10*np.log10((log_d-log_y)**2),"r", label="e - error [dB]")
plt.legend(); plt.tight_layout(); plt.show()
