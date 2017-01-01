from alex_vsd_analysis_tools import *
from pyqtgraph.Qt import QtCore, QtGui
from skimage.measure import block_reduce

app = QtGui.QApplication([])    #Starts up the QtGui; makes gl object that needs to be passed to graph routines...


#******** LOAD FILES ***********

#main_dir = '/media/cat/12TB/in_vivo/tim/alex/vsd_june_8/'
#main_dir = '/media/cat/12TB/in_vivo/tim/alex/vsd_june_24/'
#main_dir = '/media/cat/12TB/in_vivo/tim/alex/vsd_june_27/'
#main_dir = '/media/cat/12TB/in_vivo/tim/alex/vsd_june_27/Glusniff/'
main_dir = '/media/cat/12TB/in_vivo/tim/alex/vsd_june_27/VSD/'

#files=['1','2','3','4','5','6','7','8','9','10','11']

files = glob.glob(main_dir+"*.tif")

#files_picked = ['29','30','31','32']
#temp_files = []
#for file_ in files:
    #for file_picked in files_picked:
        #if file_picked in file_: temp_files.append(file_)
#files_= temp_files

#****** CHOOSE DIM REDUCTION METHOD ******
methods = ['tSNE']

files=files[int(len(files)*3/4.):]

#********* PROCESS DATA FILES **********
if False: process_data(files, main_dir, methods)

#****** LOAD PROCESSED DATA ********
if False: 
    for file_ in files[int(len(files)*2./3.):]:
        print file_
        images_filtered = np.load(file_[:-4]+'_filtered.npy')
        #DFF Computation
        data = []
        print "...computing DFF..."
        print "...cat's method..."
        images_filtered1 = images_filtered - np.nanmin(images_filtered) #+ images_raw[300]
        baseline = np.average(images_filtered1[300:], axis=0)       #Cat's method: uses minimum offset for processing
        data = (images_filtered1[300:]-baseline)/baseline
        
        print data.shape
        data = block_reduce(data, block_size=(1,2,2), func=np.mean)
        print data.shape

        #Vectorize data
        temp_data = data.reshape(data.shape[0],4096)
        print temp_data.shape

        #temp_data=temp_data[:100]
        dists = scipy.spatial.distance.pdist(temp_data,'euclidean')
        #print dists[0:10]
        
        np.save(file_[:-4]+'_'+str(len(temp_data))+'_all', dists)
        #quit()

#***** PROCESS PAIRWISE DATA - ALL DIMENSIONS ************
if False: 
    methods = ['all']
    count_pairwise_distances_all(methods, files, main_dir)

#******** PLOT 3D DATA ***************
if True: plot_pyqt_3d(app, methods, files); quit()


#******** PROCESS DIM REDUCED DATA ************
count_pairwise_distances(methods, files, main_dir)

count_velocities(methods, files, main_dir)


quit()






#count_cube_distribution(methods, files, main_dir)




