import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import scipy.ndimage as ndimage
import scipy
import os
from mpl_toolkits.mplot3d import Axes3D
from sklearn import *
import glob
import math
from math import pi

#*******************************

def dim_reduction(matrix_in, method, filename):
    
    methods = ['MDS - SMACOF', 't-SNE', 'PCA', 'Sammon']
    
    print "Computing dim reduction, size of array: ", np.array(matrix_in).shape
    
    if method==0:
        #MDS Method - SMACOF implementation Nelle Varoquaux
        #if os.path.exists(filename+'_MDS.npy')==False:
            
            array_1D = np.zeros((len(matrix_in), len(matrix_in[0].ravel())), dtype = np.float32)
            for k in range(len(matrix_in)):
                array_1D[k] = matrix_in[k].ravel()
            matrix_in = array_1D
            
            print "... MDS-SMACOF..."
            print "... pairwise dist ..."
            dists = metrics.pairwise.pairwise_distances(matrix_in)
            adist = np.array(dists)
            amax = np.amax(adist)
            adist /= amax
            
            print "... computing MDS ..."
            mds_clf = manifold.MDS(n_components=3, metric=True, n_jobs=-1, dissimilarity="precomputed", random_state=6)
            results = mds_clf.fit(adist)
            Y = results.embedding_ 

            np.save(filename+'_MDS', Y)
        #else:
            Y = np.load(filename+'_MDS.npy')
                
    elif method==1:
        ##t-Distributed Stochastic Neighbor Embedding; Laurens van der Maaten
        #if os.path.exists(filename+'_tSNE.npy')==False:
            print "... tSNE ..."
            print "... pairwise dist ..."
            
            array_1D = np.zeros((len(matrix_in), len(matrix_in[0].ravel())), dtype = np.float32)
            for k in range(len(matrix_in)):
                array_1D[k] = matrix_in[k].ravel()
            matrix_in = array_1D
            
            dists = metrics.pairwise.pairwise_distances(matrix_in)
            
            adist = np.array(dists)
            amax = np.amax(adist)
            adist /= amax
            
            print "... computing tSNE ..."
            model = manifold.TSNE(n_components=3, init='pca', random_state=0)
            Y = model.fit_transform(adist)
            #Y = model.fit(adist)
        
            np.save(filename+'_tSNE', Y)
        
        #else:
        #    Y = np.load(filename+'_tSNE.npy')

    elif method==2:

        #if os.path.exists(filename+'_PCA.npy')==False:
            print "...computing PCA..."

            array_1D = np.zeros((len(matrix_in), len(matrix_in[0].ravel())), dtype = np.float32)
            for k in range(len(matrix_in)):
                array_1D[k] = matrix_in[k].ravel()
            matrix_in = array_1D

            Y, X = PCA(matrix_in, 3)

            np.save(filename+'_PCA', Y)
        #else:
            Y = np.load(filename+'_PCA.npy')
            
                
    elif method==3:

        #if os.path.exists(filename+'_tSNE_barnes_hut.npy')==False:
            print "... computing Barnes-Hut tSNE..."
            
            array_1D = np.zeros((len(matrix_in), len(matrix_in[0].ravel())), dtype = np.float32)
            for k in range(len(matrix_in)):
                array_1D[k] = matrix_in[k].ravel()
            matrix_in = array_1D
            
            
            Y = bh_sne(np.array(matrix_in))
        
            np.save(filename+'_tSNE_barnes_hut', Y)
        #else:
            Y = np.load(filename+'_tSNE_barnes_hut.npy')

    return Y


def PCA(X, n_components):
    from sklearn import decomposition

    pca = decomposition.PCA(n_components)
    pca.fit(X)
    X=pca.transform(X)

    coords = []
    for i in range(len(X)):
         coords.append([X[i][0], X[i][1], X[i][2]])
    
    return X, np.array(coords).T #THIS IS REDUNDANT... REDUCE IT

def KMEANS(data, n_clusters):

    from sklearn import cluster, datasets
    clusters = cluster.KMeans(n_clusters, max_iter=1000, n_jobs=-1, random_state = 1032)
    clusters.fit(data)
    
    return clusters.labels_
    
    
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

def quick_mask_single_allframe(data, midline_mask_n_pixels):
    
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
    temp_array = np.ma.array(np.zeros((len(data), n_pixels,n_pixels),dtype=np.float32), mask=True)
    for k in range(len(data)):
        temp_array[k] = np.ma.masked_array(data[k], mask=generic_mask_indexes, fill_value=0)
    
    return temp_array


def PointsInCircum(r,n=100):
    return [(math.cos(2*pi/n*x)*r,math.sin(2*pi/n*x)*r) for x in xrange(0,n+1)]
          

def define_area_circle(event):
    
    global area_coords, fig, cid, circle_size
    
    if event.inaxes is not None:
        #Define N points on a circle centred at mouse click; shift circle to location
        points = np.vstack(PointsInCircum(circle_size,n=20))
        points = points + [int(event.ydata), int(event.xdata)]

        area_coords = []
        for i in range(len(points)):
            area_coords.append((points[i][0], points[i][1]))
            
        #Plot recent area
        temp_coords = np.array(area_coords)
        plt.plot(temp_coords[:,1],temp_coords[:,0],color='white',linewidth=4)
        
        fig.canvas.draw()
        
    else:
        plt.close()
        fig.canvas.mpl_disconnect(cid)
        
    return
    
def Define_cortical_areas(file_name):
    print "Defining cortical areas"

    area_names = ['limb', 'barrel','retrosplenial','visual', 'motor'] 

    sides = ['left','right']

    global area_coords, ax, fig, cid, circle_size #Not sure need all these vars
    n_pixels = 64
  
    depth = 'cortex'#, 'subcortical']
    
    #circle_sizes = [15,15,15,15,15,15,15,15]
    circle_sizes = [3,3,3,3,3,3,3,3]    #In # of pixels
    
    dir_path, filename = os.path.split(file_name) 
    
    #for depth in depths:
    counter=0
    for ctr, area in enumerate(area_names):
        circle_size = circle_sizes[ctr]
        
        for side in sides:
            save_file = dir_path+'/'+depth+"_"+area+"_"+ side
            if (os.path.exists(save_file+'.npy')==False):
                
                #images_temp = Load_max_map(file_dir, depth+ " " + area+" " + side)
                images_temp = []
                for p in range(10):
                    images_temp.append(np.load(file_name, mmap_mode='c')[0,:,5760+p*64:5760+(p+1)*64])
                images_temp = np.mean(images_temp, axis=0)
                images_temp = quick_mask_single_frame(images_temp, 2)
                #plt.imshow(img_out)
                #plt.show()

                area_coords = []
                fig, ax = plt.subplots()
                ax.imshow(images_temp)
                ax.set_title("\nDefine Location of "+depth+" " + area+" "+side, fontsize=30)
                #cid = fig.canvas.mpl_connect('button_press_event', define_area_manual)
                
                cid = fig.canvas.mpl_connect('button_press_event', define_area_circle)
                
                #fig.canvas.update()
                figManager = plt.get_current_fig_manager()
                figManager.window.showMaximized()
                plt.ylim(n_pixels,0)
                plt.xlim(0,n_pixels)
                plt.show()

                #print area_coords

                #Convert coords into x and y data for plotting; select each column, then tack on first coordinate to complete circle
                area_coords.append(area_coords[0])
                area_coords = np.array(area_coords)
                
                from matplotlib import path
                #Compute pixel locations inside cropped area
                p = path.Path(area_coords)
                all_pts = []
                for i in range(64):
                    for j in range(64):
                        all_pts.append([i,j])
                pts_inside = p.contains_points(all_pts)
                
                #Generate mask for saving 
                mask_save = np.zeros((64,64),dtype=np.int8)+1
                for i in range(64):
                    for j in range(64):
                        if pts_inside[i*64+j]==True:
                            mask_save[i,j]=False
                
                #Save mask
                np.save(save_file, mask_save)
                #np.save(file_name + "subcortical_"+area+"_"+ side, mask_save)
                
                ##Save contour
                #np.save(save_file+'_contour', area_coords)
                #np.save(file_dir + "subcortical_"+area+"_"+ side+'_contour', mask_save, area_coords)
                
        counter+=1

def plot_data(data, colours, original_image, methods, method, n_clusters, labels, file_name, mid_line_mask, block, unit):


    print file_name
    path, file_name = os.path.split(file_name) 
    if True:
        #***************************************************************************
        #PLot rasters
        ax = plt.subplot(111)
        
        raster_file = glob.glob(path+"/"+"unit_"+str(unit).zfill(2)+"*imagingspikes.txt")[0]
        spikes = np.loadtxt(raster_file)
        spikes = spikes - np.min(spikes)

        ymin = np.zeros(len(labels),dtype=np.float32)-labels*5
        ymax = np.zeros(len(labels),dtype=np.float32)-(labels+1)*5+1
        print len(spikes), len(labels)
        
        spikes = spikes[:len(labels)]   #Trim the last few spikes that weren't computed in the STM routine due to parallelization round off
        colors_out = []
        for k in range(len(labels)):
            colors_out.append(colours[labels[k]])
        plt.vlines(spikes, ymin, ymax, color=colors_out, linewidth=1, alpha=0.5)
        
        #Plot original rasters
        offset=np.max(labels)+1
        ax.yaxis.set_ticks([])

        ymin = np.zeros(len(labels),dtype=np.float32)-offset*5
        ymax = np.zeros(len(labels),dtype=np.float32)-(offset+1)*5+1
        plt.vlines(spikes, ymin, ymax, color='black', linewidth=1, alpha=0.5)

        plt.xlabel("Time (sec)", fontsize = 30)
        plt.ylabel("All Spikes  Sub Groups", fontsize = 30)
        plt.tick_params(axis='both', which='both', labelsize = 30)
        plt.xlim(0, np.max(spikes)+5)
        plt.show()

        #*****************************************************************************
        #Plot ISI histograms
        for k in range(n_clusters):
            #ax = plt.subplot(1,n_clusters,k+1)
            ax = plt.subplot(2,2,k+1)
            ax.yaxis.set_ticks([])
            indexes = np.where(labels==k)[0]
            temp_raster = spikes[indexes]
            isi = temp_raster[1:]-temp_raster[:-1]
            
            bin_width = 0.005
            y = np.histogram(isi, bins = np.arange(0,0.2,bin_width))
            plt.bar(y[1][:-1], y[0], bin_width, color=colours[k])        
            plt.ylim(0,np.max(y[0]))
            
            plt.xlabel("Time (sec)", fontsize = 30)
            plt.tick_params(axis='both', which='both', labelsize=20)
        
        plt.suptitle("ISI Histograms", fontsize = 30)
        plt.show()
        
        #PLot all spikes ISI histogram
        ax = plt.subplot(2,2,1)
        isi = spikes[1:]-spikes[:-1]
        bin_width = 0.005
        y = np.histogram(isi, bins = np.arange(0,0.2,bin_width))
        plt.bar(y[1][:-1], y[0], bin_width, color='black')        
        plt.ylim(0,np.max(y[0]))
        plt.show()
        
        #*************************************************************************************
        #Plot partitioned rasters
        ymin = np.zeros(len(labels),dtype=np.float32)-labels*5
        ymax = np.zeros(len(labels),dtype=np.float32)-(labels+1)*5+1
        print len(spikes), len(labels)
        
        spikes = spikes[:len(labels)]   #Trim the last few spikes that weren't computed in the STM routine due to parallelization round off
        colors_out = []
        for k in range(len(labels)):
            colors_out.append(colours[labels[k]])
        plt.vlines(spikes, ymin, ymax, color=colors_out, linewidth=1, alpha=0.5)
        
        #Plot original rasters
        offset=np.max(labels)+1
        ax.yaxis.set_ticks([])

        ymin = np.zeros(len(labels),dtype=np.float32)-offset*5
        ymax = np.zeros(len(labels),dtype=np.float32)-(offset+1)*5+1
        plt.vlines(spikes, ymin, ymax, color='black', linewidth=1, alpha=0.5)

        plt.xlabel("Time (sec)", fontsize = 30)
        plt.ylabel("All Spikes  Sub Groups", fontsize =30)
        plt.tick_params(axis='both', which='both', labelsize=30)
        plt.xlim(0, np.max(spikes)+5)
        plt.suptitle("Partitioned Rasters", fontsize=30)
        plt.show()
        
        
        
        #******************************************************************************
        #Plot Scatter Distributions
        if True: 
            fig = plt.figure()
            ax = fig.add_subplot(111)

            for k in range(len(data)): 
                ax.scatter(data[k][0], data[k][1], c=colours[labels[k]], marker='o')

            plt.show()
            

        if False: 
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            for k in range(len(data)): 
                ax.scatter(data[k][0], data[k][1], data[k][2], c=colours[labels[k]], marker='o')

            ax.set_xlabel('X Label')
            ax.set_ylabel('Y Label')
            ax.set_zlabel('Z Label')

            plt.show()
            
        
        #******************************************************************************
        #Plot motifs
        print "... plotting ..."
        for k in range(n_clusters):
            print "...cluster: ", k
            ax=plt.subplot(n_clusters+1,1,k+1)
            indexes = np.where(labels==k)[0]
            temp_ave = np.mean(original_image[indexes], axis=0)
            
            img_stack = []
            for p in range(temp_ave.shape[1]/64):
                img_stack.append(temp_ave[:,p*64:(p+1)*64])
            
            temp_stack = []
            for p in range(0, len(img_stack), block):
                temp_stack.append(np.mean(img_stack[p:p+block], axis=0))

            temp_stack = quick_mask_single_allframe(temp_stack, mid_line_mask)
            
            img_out = np.ma.hstack((temp_stack))
            img_out[:, img_out.shape[1]/2-2:img_out.shape[1]/2+2] = np.min(img_out)
            
            v_max = np.nanmax(np.abs(img_out)); v_min = -v_max

            plt.imshow(img_out, vmin=v_min, vmax=v_max)
            
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            
            ax.set_title("Cluster: "+colours[k]+" # Spks: "+str(len(indexes)) + " DF/F (max/min): " + str(round(v_max*100,1)) +"%", fontsize =30)

        #Add total average motif
        ax=plt.subplot(n_clusters+1, 1, n_clusters+1)
        indexes = np.where(labels==k)[0]
        temp_ave = np.mean(original_image, axis=0)
        
        img_stack = []
        for p in range(temp_ave.shape[1]/64):
            img_stack.append(temp_ave[:,p*64:(p+1)*64])
        
        temp_stack = []
        for p in range(0, len(img_stack), block):
            temp_stack.append(np.mean(img_stack[p:p+block], axis=0))

        temp_stack = quick_mask_single_allframe(temp_stack, mid_line_mask)
        
        img_out = np.ma.hstack((temp_stack))
        v_max = np.nanmax(np.abs(img_out)); v_min = -v_max

        img_out[:, img_out.shape[1]/2-2:img_out.shape[1]/2+2] = -100


        plt.imshow(img_out, vmin=v_min, vmax=v_max)
        
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        plt.title("Average of All Spikes: "+str(original_image.shape[0]) + "   DF/F (max/min): " + str(round(v_max*100,1)) +"%", fontsize = 30)

        plt.suptitle(file_name + "\nDimension Reduction Method: " + methods[method], fontsize = 30)
        plt.show()    
    
    #****************************************************************************
    #Track max val location at specific ROI
    depth = 'cortex'
    areas = ['motor', 'retrosplenial', 'barrel']
    side = 'left'

    if False: 
        for area in areas: 
            mask_file = path + '/' + depth + "_" + area + "_" + side       #Load ROI Mask
            temp_data = np.load(mask_file+'.npy')
            mask_stack_ravel = temp_data.ravel()
            print mask_stack_ravel.shape

            mask_indexes = np.where(mask_stack_ravel==0)[0]      #Pick only values not masked, i.e. mask = False
            
            #********************
            #Plot motifs
            overall_max = []
            print "... plotting ..."
            for k in range(n_clusters):
                print "...cluster: ", k
                #ax=plt.subplot(n_clusters+1,1,k+1)
                ax=plt.subplot(1,1,1)
                
                #Select clustered motifs only
                indexes = np.where(labels==k)[0]
                temp_ave = np.mean(original_image[indexes], axis=0)

                #Convert back to multi-frame stack
                img_stack = []
                for p in range(temp_ave.shape[1]/64):
                    img_stack.append(temp_ave[:,p*64:(p+1)*64])
                
                img_stack = np.array(img_stack)
                print img_stack.shape

                #Find location of max value - within ROI; NB: THIS MAY NOT ALWAYS WORK!!!
                max_loc = []; max_val = []
                for p in range(len(img_stack)):
                    temp_max = np.argmax(np.abs(img_stack[p].ravel()[mask_indexes])) #Find the abs_max value in each ROI chosen
                    max_loc.append(temp_max)
                    max_val.append(img_stack[p].ravel()[mask_indexes][temp_max])
                
                max_val_argument = np.argmax(np.abs(max_val))
                max_val = max_val[max_val_argument]*1E2
                max_loc = max_loc[max_val_argument]

                print "...max locations: ", max_loc, " max_val: ", max_val
                overall_max.append(max_val)

                #Track max/min value across all frames
                trace_max = []
                for p in range(img_stack.shape[0]):                     #Loop over all frames
                    trace_max.append(img_stack[p].ravel()[mask_indexes][max_loc])
                trace_max = np.array(trace_max)*1E2

                plt.plot(trace_max, color=colours[k], linewidth=3)
                
                np.savetxt(path+'/'+str(unit).zfill(2)+"_"+area+"_cluster"+str(k), trace_max)
            
            #********************
            #Plot motifs
            print "... plotting ..."
            
            #Select clustered motifs only
            temp_ave = np.mean(original_image, axis=0)

            #Convert back to multi-frame stack
            img_stack = []
            for p in range(temp_ave.shape[1]/64):
                img_stack.append(temp_ave[:,p*64:(p+1)*64])
            
            img_stack = np.array(img_stack)
            print img_stack.shape

            #Find location of max value - within ROI CHOSEN ABOVE
            max_loc = []; max_val = []
            for p in range(len(img_stack)):
                temp_max = np.argmax(np.abs(img_stack[p].ravel()[mask_indexes])) #Find the abs_max value in each ROI chosen
                max_loc.append(temp_max)
                max_val.append(img_stack[p].ravel()[mask_indexes][temp_max])
            
            max_val_argument = np.argmax(np.abs(max_val))
            max_val = max_val[max_val_argument]*1E2
            max_loc = max_loc[max_val_argument]
            
            print "...max locations: ", max_loc, " max_val: ", max_val

            #Track max/min value across all frames
            trace_max = []
            for p in range(img_stack.shape[0]):                     #Loop over all frames
                trace_max.append(img_stack[p].ravel()[mask_indexes][max_loc])
            trace_max = np.array(trace_max)*1E2
                    
            plt.plot(trace_max, color='black', linewidth=3)
                
            np.savetxt(path+'/'+str(unit).zfill(2)+"_"+area+"_allspikes"+str(k), trace_max)


            plt.plot([0,180],[0,0], 'r--', color='black', linewidth=3, alpha=0.5)
            plt.plot([90,90],[-10,10], 'r--', color='black', linewidth=3, alpha=0.5)
            
            plt.ylim(int(-np.max(np.abs(overall_max))-1),int(np.max(np.abs(overall_max))+1))
            plt.xlim(0,180)

            plt.ylabel("DF/F (%)", fontsize=30)
            
            plt.title("STMTD for Left-"+area+" Cortex", fontsize = 30)
            plt.xlabel("Time (sec)", fontsize=30)

            old_xlabel = np.arange(0, len(trace_max)+3,30)
            new_xlabel = np.arange(-3.0,3.1, 1)
            plt.xticks(old_xlabel, new_xlabel, fontsize=30)    

            plt.tick_params(axis='both', which='both', labelsize=30)


            plt.show()
                

    #********************************************************************
    #Plot 3D trace DF/F distributions
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    all_data = True

    #Plot partitioned dynamics
    for k in range(4):
        temp_traces = []
        for area in areas: 
            temp_traces.append(np.loadtxt(path+'/'+str(unit).zfill(2)+"_"+area+"_cluster"+str(k)))

        temp_traces = np.vstack(temp_traces)
        
        if all_data:    #Plot -1..+1sec
            ax.plot(temp_traces[0][len(temp_traces[0])/2-30:len(temp_traces[0])/2+30], 
                    temp_traces[1][len(temp_traces[0])/2-30:len(temp_traces[0])/2+30], 
                    temp_traces[2][len(temp_traces[0])/2-30:len(temp_traces[0])/2+30], c=colours[k], marker='o')
        else: #Plot everything
            ax.plot(temp_traces[0], 
                    temp_traces[1], 
                    temp_traces[2], c=colours[k], marker='o')
                

    #Plot all spike dynamics
    temp_traces = []
    for area in areas: 
        temp_traces.append(np.loadtxt(path+'/'+str(unit).zfill(2)+"_"+area+"_allspikes"+str(k)))

    temp_traces = np.vstack(temp_traces)
    
    if all_data:    #Plot -1..+1sec
        ax.plot(temp_traces[0][len(temp_traces[0])/2-30:len(temp_traces[0])/2+30], 
                temp_traces[1][len(temp_traces[0])/2-30:len(temp_traces[0])/2+30], 
                temp_traces[2][len(temp_traces[0])/2-30:len(temp_traces[0])/2+30], c='black', marker='o', linewidth=5)
    else: #Plot everything
        ax.plot(temp_traces[0], 
                temp_traces[1], 
                temp_traces[2], c='black', marker='o', linewidth=5)
                



    ax.set_xlabel("DF/F (%): "+areas[0], labelpad=30, fontsize=30)
    ax.set_ylabel("DF/F (%): "+areas[1], labelpad=30, fontsize=30)
    ax.set_zlabel("DF/F (%): "+areas[2], labelpad=30, fontsize=30)

    plt.tick_params(axis='both', which='both', labelsize=20)
    #ax.tick_params(axis='x', pad=-10)
    #ax.tick_params(axis='y', pad=-10)
    #ax.tick_params(axis='z', pad=-10)


    plt.show()
            


def partition_data_roi(base_filename, cells, methods, methods_array, colours, area_names, sides, mid_line_mask, block, crop_method, unit, n_clusters):
    
    c = unit
    print "... processing cell: ", c

    filename = glob.glob(base_filename+str(c).zfill(2)+"*spikes.npy")
    if len(filename)!=1: 
        print "...skipping cell: "
        return

    filename = filename[0]
    
    img_stack = np.load(filename, mmap_mode='c')
    print img_stack.shape
        
    #************************************************************************    
    #Plot example motifs
    motifs = [1]
    motifs.append(int(len(img_stack)/5)*2+25)
    motifs.append(int(len(img_stack)/5)*3-1)
    motifs.append(int(len(img_stack)/5)*4-1)
    motifs.append(int(len(img_stack)/5)*5-1)

    fig = plt.figure()

    if True:
        for k in range(len(motifs)):
            print "...motif: ", motifs[k]
            ax=plt.subplot(5,1,k+1)
            temp_ave = img_stack[motifs[k]]
            
            temp_stack = []
            for p in range(temp_ave.shape[1]/64):
                temp_stack.append(temp_ave[:,p*64:(p+1)*64])
            
            out_stack = []
            for p in range(0, len(temp_stack), block):
                out_stack.append(np.mean(temp_stack[p:p+block], axis=0))

            temp_stack = quick_mask_single_allframe(out_stack, mid_line_mask)
            
            img_out = np.ma.hstack((temp_stack))
            img_out[:, img_out.shape[1]/2-2:img_out.shape[1]/2+2] = np.min(img_out)
            
            v_max = np.nanmax(np.abs(img_out)); v_min = -v_max

            plt.imshow(img_out, vmin=v_min, vmax=v_max)
            
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            
            ax.set_title("Spike #: "+str(motifs[k]) + " DF/F (max/min): " + str(round(v_max*100,1)) +"%", fontsize =30)
    plt.suptitle(filename+"\nExample Single Spike Motifs", fontsize = 30)
    plt.show()

    global area_coords, ax, fig, cid, circle_size #Not sure need all these vars
    n_pixels = 64
  
    depth = 'cortex'#, 'subcortical']
    
    #for depth in depths:
    counter=0
    
    dir_path, file_name = os.path.split(filename) 

    vector_stack = []
    for ctr, area in enumerate(area_names):
        #circle_size = circle_sizes[ctr]
        
        for side in sides:
            
            save_file = filename[:-4] + '_' + depth + "_" + area + "_" + side+"_roi"
            out_file = save_file+"_"+crop_method

            if os.path.exists(out_file+'.npy')==False:
                
                mask_file = dir_path + '/' + depth + "_" + area + "_" + side       #Load ROI Mask
                
                temp_data = np.load(mask_file+'.npy')
                mask_stack_ravel = temp_data.ravel()
                print mask_stack_ravel.shape

                indexes = np.where(mask_stack_ravel==0)[0]      #Pick only values not masked, i.e. mask = False
                
                area_vector_stack= []
                print "... area: ", area, " side: ", side
                for s in range(img_stack.shape[0]):                     #Loop over all spikes
                    area_vector = []
                    for k in range(img_stack.shape[2]/64):              #Loop over all frames for each spike

                        img_frame = img_stack[s, :, k*64:(k+1)*64]

                        img_frame_ravel = img_frame.ravel()
                        
                        if crop_method == 'ave': area_vector.append(np.mean(img_frame_ravel[indexes]))
                        if crop_method == 'max': area_vector.append(np.max(img_frame_ravel[indexes]))

                    area_vector = np.hstack(area_vector)
                        
                    #Normalize the spike stack
                    if True: 
                        val_max = np.nanmax(area_vector); val_min = np.nanmin(area_vector)
                        area_vector = (area_vector - val_min)/(val_max-val_min)

                    area_vector_stack.append(area_vector)
                
                    #plt.plot(np.hstack(area_vector))
                area_vector_stack = np.vstack(area_vector_stack)
                
                #print area_vector_stack.shape

                np.save(out_file, area_vector_stack)
                vector_stack.append(area_vector_stack)
    
            else:
                vector_stack.append(np.load(out_file+'.npy'))
    
    vector_stack = np.array(vector_stack)
    vector_stack = np.swapaxes(vector_stack, 0, 1)
    #vector_stack = np.swapaxes(vector_stack, 1, 2)
    print vector_stack.shape
   
    
    all_stack = []
    for p in range(len(vector_stack)):
        all_stack.append([])
        for r in range(len(vector_stack[p])):
            #all_stack[p].extend(vector_stack[p][r])
            #all_stack[p].extend(vector_stack[p][r][vector_stack.shape[2]/2:])
            all_stack[p].extend(vector_stack[p][r][vector_stack.shape[2]/2:vector_stack.shape[2]/2+10])
            
            #temp_mean = np.average(vector_stack[p][r][vector_stack.shape[2]/2:vector_stack.shape[2]/2+6], axis=0)
            #all_stack[p].extend([temp_mean])
    
    data_out = np.array(all_stack)
        
           
    #ROI Dimension Reduction
    for method in methods_array: 
        path_name, file_name = os.path.split(filename)
         
        #data = dim_reduction(original_image, method, filename)

        data = dim_reduction(data_out, method, filename+"_roi_dimreduction")

        #n_clusters = 4

        #Clustering Step
        labels = KMEANS(data, n_clusters)
        
        #Plot Traces Being Clustered
        if False: 
            ave_traces = []
            for k in range(4): ave_traces.append([])
            for k in range(len(data_out)):
                #plt.plot(np.hstack(np.swapaxes(vector_stack[k], 0, 1)), color=colours[labels[k]], alpha=0.3)
                plt.plot(data_out[k], color = colours[labels[k]], alpha =0.1)
                ave_traces[labels[k]].append(data_out[k])
            
            for k in range(4):
                plt.plot(np.mean(ave_traces[k], axis=0), color=colours[k], linewidth=3)
            plt.show()
            
        #Test by randomizing labels
        if False:
            np.random.shuffle(labels)

        if True: 
            plot_data(data, colours, img_stack, methods, method, n_clusters, labels, filename, mid_line_mask, block, unit)
            
             
           
def partition_data(base_filename, cells, methods, methods_array, colors):
    
    for c in cells: 
        print "... processing cell: ", c

        filename = glob.glob(base_filename+str(c).zfill(2)+"*spikes.npy")
        if len(filename)!=1: 
            print "...skipping cell: "
            continue

        filename = filename[0]
        for method in methods_array: 
            path_name, file_name = os.path.split(filename)
             
            original_image = np.load(filename, mmap_mode='c')

            if original_image.shape[0] < 250: continue
            
            print "... # spikes: ", original_image.shape[0]

            print "... dim reduction method: ", methods[method]

            #array_1D = np.zeros(original_image.shape, dtype = np.float32)
            #array_1D = np.zeros((len(original_image), 729088), dtype = np.float32)
            #array_1D = np.zeros((len(original_image), 737280), dtype = np.float32)
            #for k in range(len(original_image)):
            #    array_1D[k] = original_image[k].ravel()

            #data = dim_reduction(original_image, method, filename)
            data = dim_reduction(original_image[:,:, original_image.shape[2]/2:], method, filename+"_half_stack")
            #data = dim_reduction(original_image[:,:, original_image.shape[2]/2:original_image.shape[2]/2+64*30], method, filename+"_1sec_stack")

            n_clusters = 4

            labels = KMEANS(data, n_clusters)

            if plotting: 
                plot_data(data, colours, original_image, methods, method, mid_line_mask)
                
                
                
                
#********************************

colours = ['green','red','blue', 'magenta','purple','orange','cyan','pink','grey', 'orange']

plotting = True

#method = 1

#base_filename = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-14-allelectrodeinthalamus-is0/stack1D_2015-12-2-14-allelectrodeinthalamus-is0_unit'
#base_filename = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-2/stack1D_2015-7-22-2_unit'
#base_filename = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-5/stack1D_2015-7-22-5_unit'
#base_filename = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-6-10electrodeincortex-iso0/stack1D_2015-12-2-6-10electrodeincortex-iso0_unit'
#base_filename = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-18/2015-11-18-8-9electrodein-iso0/stack1D_2015-11-18-8-9electrodein-iso0_unit'
#base_filename = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-5-10electrodein-iso0/stack1D_2015-11-27-5-10electrodein-iso0_unit'
base_filename = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-15-allelectrodeinthalamus-is0/stack1D_2015-12-2-15-allelectrodeinthalamus-is0_unit'
#base_filename = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-12-allelectrodeinthalamus-iso1/stack1D_2015-12-2-12-allelectrodeinthalamus-iso1_unit'
#base_filename = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-11/2015-12-11-12-allelectinthalamus-iso1/stack1D_2015-12-11-12-allelectinthalamus-iso1_unit'

#filename = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-6-10electrodeincortex-iso0/stack1D_2015-12-2-6-10electrodeincortex-iso0_unit02_ch09_all_3sec_window_00283_spikes.npy'
#filename = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-6-10electrodeincortex-iso0/stack1D_2015-12-2-6-10electrodeincortex-iso0_unit05_ch13_all_3sec_window_01936_spikes.npy'
#filename = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-6-10electrodeincortex-iso0/stack1D_2015-12-2-6-10electrodeincortex-iso0_unit09_ch14_all_3sec_window_01728_spikes.npy'
#filename = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-6-10electrodeincortex-iso0/stack1D_2015-12-2-6-10electrodeincortex-iso0_unit08_ch14_all_3sec_window_01886_spikes.npy'

   
#area_names = ['limb', 'barrel','retrosplenial','visual', 'motor'] 
area_names = ['barrel', 'retrosplenial', 'motor'] 
sides = ['left']#,'right']
crop_method = 'ave'  #,['max','ave']


mid_line_mask = 2
block = 10 #30Hz recs
#block = 15 #50Hz recs

n_clusters = 4
        
cells = np.arange(0,100,1)
#unit = 1 #7-22-2
#unit = 16 #7-22-5
#unit = 4 #12-2; cells: 7, 12, 8, 9, 13, 4, 5
#unit = 2 #11-8-9
#unit = 12 #11-27-5
unit = 33 #2015-12-215
#unit = 65 #12-2-12
#unit = 2 #12-11-12

methods = ['MDS - SMACOF', 't-SNE', 'PCA', 'tSNE-Barnes_hut']
methods_array = np.arange(0,3,1)
methods_array = [2]


temp_name = glob.glob(base_filename+str(unit).zfill(2)+"*spikes.npy")
filename = temp_name[0]


Define_cortical_areas(filename)


#partition_data(base_filename, cells, methods, methods_array, colors)


partition_data_roi(base_filename, cells, methods, methods_array, colours, area_names, sides, mid_line_mask, block, crop_method, unit, n_clusters)




quit()




#Test normality of pixel distributions
filename = glob.glob(base_filename+str(cells[0]).zfill(2)+"*spikes.npy")[0]

original_image = np.load(filename, mmap_mode='c')

print original_image.shape

from scipy.stats import kstest
import statsmodels.api as sm

vertical_pixels = np.arange(0,64,1)

frames = np.arange(-10,30,1)
for frame in frames:
    print "...frame: ", frame
    
    start_frame_pixel = 5760+64*frame
    horizontal_pixels = np.arange(start_frame_pixel,start_frame_pixel+64,1)

    out_img = []
    for v in vertical_pixels:
        h_lines = []
        
        for h in horizontal_pixels:
            data = original_image[:,v,h]
            
            #bin_width = 0.005   # histogram bin width in usec
            #y = np.histogram(data, bins = np.arange(-.2,.2,bin_width))
            ##plt.bar(y[1][:-1]*1E-3, y[0])/np.max(y[0])*len(locked_spikes), bin_width*1E-3, color='blue', alpha=0.2)
            #plt.bar(y[1][:-1], y[0], bin_width, color='blue')

            #plt.title(str(scipy.stats.mstats.normaltest(data, axis=0)[1]))
            #plt.title(str(kstest(data,'norm')))

            W, p = scipy.stats.shapiro(data)
            #h_lines.append((p + sm.stats.normal_ad(data)[1])/2.)
            h_lines.append(p)

        out_img.append(np.hstack(h_lines))

    out_img = np.vstack(out_img)

    out_img = quick_mask_single_frame(out_img, 5)
    
    ax = plt.subplot(6,7,frame+11)
        
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

    plt.imshow(out_img,vmin = 0, vmax=0.01, interpolation='none')

plt.show()

















