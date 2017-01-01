
#Mouse class and subclasses for processing lever pull data projec Greg Silasi and Tim Murphy lab
#Code: Cat Mitelut

import os
import glob
import numpy as np
import struct
import string, re
import scipy
import tifffile as tiff
from skimage.measure import block_reduce
import matplotlib.pyplot as plt
import shutil
from imreg import *
#from lever_pull_utils import *
from numpy import inf


class Ptcs(object):
    """Polytrode clustered spikes file neuron record"""
    def __init__(self, file_name):
        
        f = open(file_name, "rb")
        self.sorted_file = file_name
        self.name = file_name
        self.full_path =file_name
        # call the appropriate method:
        self.VER2FUNC = {1: self.readHeader, 2: self.readHeader, 3: self.readHeader}

        self.readHeader(f)
        
        self.nid = []  #Make unique unit id list for loading later.
        
        self.loadData(self.nsamplebytes, f)
        
        f.close()

    def __getstate__(self):
        """Instance methods must be excluded when pickling"""
        d = self.__dict__.copy()
        try: del d['VER2FUNC']
        except KeyError: pass
        return d

    def readHeader(self, f):
        """Read in neuron record of .ptcs file version 3. 'zpos' field was replaced
        by 'sigma' field.
        nid: int64 (signed neuron id, could be -ve, could be non-contiguous with previous)
        ndescrbytes: uint64 (nbytes, keep as multiple of 8 for nice alignment, defaults to 0)
        descr: ndescrbytes of ASCII text
        (padded with null bytes if needed for 8 byte alignment)
        clusterscore: float64
        xpos: float64 (um)
        ypos: float64 (um)
        sigma: float64 (um) (Gaussian spatial sigma)
        nchans: uint64 (num chans in template waveforms)
        chanids: nchans * uint64 (0 based IDs of channels in template waveforms)
        maxchanid: uint64 (0 based ID of max channel in template waveforms)
        nt: uint64 (num timepoints per template waveform channel)
        nwavedatabytes: uint64 (nbytes, keep as multiple of 8 for nice alignment)
        wavedata: nwavedatabytes of nsamplebytes sized floats
        (template waveform data, laid out as nchans * nt, in uV,
        padded with null bytes if needed for 8 byte alignment)
        nwavestdbytes: uint64 (nbytes, keep as multiple of 8 for nice alignment)
        wavestd: nwavestdbytes of nsamplebytes sized floats
        (template waveform standard deviation, laid out as nchans * nt, in uV,
        padded with null bytes if needed for 8 byte alignment)
        nspikes: uint64 (number of spikes in this neuron)
        spike timestamps: nspikes * uint64 (us, should be sorted)
        """

        self.nid = int(np.fromfile(f, dtype=np.int64, count=1)) # nid
        self.ndescrbytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # ndescrbytes
        self.descr = f.read(self.ndescrbytes).rstrip('\0 ') # descr

        if self.descr:
            try:
                self.descr = eval(self.descr) # might be a dict
            except: pass

        self.nneurons = int(np.fromfile(f, dtype=np.uint64, count=1)) # nneurons
        self.nspikes = int(np.fromfile(f, dtype=np.uint64, count=1)) # nspikes
        self.nsamplebytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # nsamplebytes
        self.samplerate = int(np.fromfile(f, dtype=np.uint64, count=1)) # samplerate
        self.npttypebytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # npttypebytes

        self.pttype = f.read(self.npttypebytes).rstrip('\0 ') # pttype

        self.nptchans = int(np.fromfile(f, dtype=np.uint64, count=1)) # nptchans
        self.chanpos = np.fromfile(f, dtype=np.float64, count=self.nptchans*2) # chanpos
        self.chanpos.shape = self.nptchans, 2 # reshape into rows of (x, y) coords
        self.nsrcfnamebytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # nsrcfnamebytes
        self.srcfname = f.read(self.nsrcfnamebytes).rstrip('\0 ') # srcfname
        # maybe convert this to a proper Python datetime object in the Neuron:
        self.datetime = float(np.fromfile(f, dtype=np.float64, count=1)) # datetime (days)
        self.ndatetimestrbytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # ndatetimestrbytes
        self.datetimestr = f.read(self.ndatetimestrbytes).rstrip('\0 ') # datetimestr


    def loadData(self, n_bytes, f):
        #call the appropriate method:
        #self.VER2FUNC = {1: self.read_ver_1, 2:self.read_ver_2, 3:self.read_ver_3}
        self.nsamplebytes = n_bytes
        self.wavedtype = {2: np.float16, 4: np.float32, 8: np.float64}[self.nsamplebytes]

        self.n_units=self.nneurons
        self.units=[None]*self.n_units
        self.uid = [None]*self.n_units  #Unique id for full track sorts
        self.n_sorted_spikes = [None]*self.n_units
        self.ptp=np.zeros((self.n_units), dtype=np.float32)
        self.size = []
        self.maxchan = []
        self.sigma = []
        self.xpos = []
        self.ypos = []
        self.wavedata = []

        for k in range(self.n_units):
            self.readUnit(f)
            self.units[k]= self.spikes

            if 'martin' in self.full_path:
                self.uid[k]= self.nid
            else: #All other sorts are from Nick's SS so should be the same
                self.uid[k]= self.nid-1
               
            #print "SAMPLERATE: ", self.samplerate
            #if ptcs_flag: #Martin's data has wrong flag for saves
            
            #CONVERT UNITS TO TIMESTEPS
            #self.units[k]=[x*self.samplerate/1E+6 for x in self.units[k]] #Converts spiketimes from usec to timesteps
            
            
            #else:
            #    self.units[k]=[x*self.samplerate/2/1E+6 for x in self.units[k]] #Converts spiketimes from usec to timesteps

            self.n_sorted_spikes[k] = len(self.units[k])
            self.size.append(self.nspikes)
            self.maxchan.append(self.maxchanu)
            self.sigma.append(self.zps)
            self.xpos.append(self.xps)
            self.ypos.append(self.yps)
            self.wavedata.append(self.wvdata)

            #self.ptp[k]=max(self.wavedata[np.where(self.chans==self.maxchanu)[0][0]]) - \
            #            min(self.wavedata[np.where(self.chans==self.maxchanu)[0][0]]) #compute PTP of template;


        f.close()


    def readUnit(self,f):
        self.nid = int(np.fromfile(f, dtype=np.int64, count=1)) # nid
        self.ndescrbytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # ndescrbytes
        self.descr = f.read(self.ndescrbytes).rstrip('\0 ') # descr

        if self.descr:
            try:
                self.descr = eval(self.descr) # might be a dict
            except: pass

        self.clusterscore = float(np.fromfile(f, dtype=np.float64, count=1)) # clusterscore
        self.xps = float(np.fromfile(f, dtype=np.float64, count=1)) # xpos (um)
        self.yps = float(np.fromfile(f, dtype=np.float64, count=1)) # ypos (um)
        self.zps = float(np.fromfile(f, dtype=np.float64, count=1)) # zpos (um) #Replaced by spatial sigma
        self.nchans = int(np.fromfile(f, dtype=np.uint64, count=1)) # nchans
        self.chans = np.fromfile(f, dtype=np.uint64, count=self.nchans) #NB: Some errors here from older .ptcs formats
        self.maxchanu = int(np.fromfile(f, dtype=np.uint64, count=1)) # maxchanid

        self.nt = int(np.fromfile(f, dtype=np.uint64, count=1)) # nt: number of time points in template

        self.nwavedatabytes, self.wvdata = self.read_wave(f) #TEMPLATE

        self.nwavestdbytes, self.wavestd = self.read_wave(f) #STANDARD DEVIATION
        self.nspikes = int(np.fromfile(f, dtype=np.uint64, count=1)) # nspikes

        self.spikes = np.fromfile(f, dtype=np.uint64, count=self.nspikes) # spike timestamps (us):

        # convert from unsigned to signed int for calculating intervals:
        #self.spikes = np.asarray(self.spikes, dtype=np.float64)

            
    def read_wave(self, f):
        """Read wavedata/wavestd bytes"""
        # nwavedata/nwavestd bytes, padded:
        nbytes = int(np.fromfile(f, dtype=np.uint64, count=1))
        fp = f.tell()
        count = nbytes // self.nsamplebytes # trunc to ignore any pad bytes
        X = np.fromfile(f, dtype=self.wavedtype, count=count) # wavedata/wavestd (uV)
        if nbytes != 0:
            X.shape = self.nchans, self.nt # reshape
        f.seek(fp + nbytes) # skip any pad bytes
        return nbytes, X

    def rstrip(s, strip):
        """What I think str.rstrip should really do"""
        if s.endswith(strip):
            return s[:-len(strip)] # strip it
        else:
            return s

    def read(self):
        self.nid = self.parse_id()
        with open(self.fname, 'rb') as f:
            self.spikes = np.fromfile(f, dtype=np.int64) # spike timestamps (us)
        self.nspikes = len(self.spikes)

    
def quick_mask_single_frame(data, midline_mask_n_pixels):
    
    n_pixels = 256
        
    generic_mask_file = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-11/genericmask.txt'
    if (os.path.exists(generic_mask_file)==True):
        generic_coords = np.int32(np.loadtxt(generic_mask_file))
    else:
        print "...generic mask not found..."
        return
    
    #Load generic mask
    generic_mask_indexes=np.zeros((n_pixels,n_pixels))
    for i in range(len(generic_coords)):
        generic_mask_indexes[int(generic_coords[i][0])][int(generic_coords[i][1])] = True

    #Load midline mask
    for i in range(midline_mask_n_pixels):
        generic_mask_indexes[:,n_pixels/2+int(midline_mask_n_pixels/2)-i]=True

    #Apply full mask; probably FASTER METHOD
    temp_array = np.ma.array(np.zeros((n_pixels,n_pixels),dtype=np.float32), mask=True)
    temp_array = np.ma.masked_array(data, mask=generic_mask_indexes, fill_value=0)
    
    return temp_array
    

#filename = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-14-allelectrodeinthalamus-is0_mua/2015-12-2-14-allelectrodeinthalamus-is0_mua_hp.ptcs'
#data = Ptcs(filename)
    
#print len(data.units)
#for k in range(len(data.units)):
    #file_out = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-14-allelectrodeinthalamus-is0_mua/'
    #file_out = file_out + "unit_"+str(k).zfill(2)+"_channel_"+str(k+1).zfill(2)+"_ptp_050.csv"
    #np.savetxt(file_out, data.units[k]*1E-6)
    

#quit()


filename = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-12-2/2015-12-2-14-allelectrodeinthalamus-is0_mua/img_avg_2015-12-2-14-allelectrodeinthalamus-is0_mua_unit15_ch16_all_3sec_window_02184_spikes.npy'

data = np.load(filename)
temp_data = []
mid_mask = 5
for k in range(len(data)):
    temp_data.append(quick_mask_single_frame(data[k], mid_mask))

v_max = np.nanmax(np.ma.abs(temp_data)); v_min = -v_max

for k in range(len(data)):
#for k in range(10):
    print "... img: ", k
    ax=plt.subplot(10,18,k+1)
    
    print v_min, v_max
    plt.imshow(temp_data[k], vmin=v_min, vmax=v_max)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
plt.show()
#print np.load('/media/cat/12TB/in_vivo/tim/yuki/AQ2/tif_files/AQ2am_Dec21_30Hz/AQ2am_Dec21_30Hz_img_rate.npy')
#quit()


#file_name = '/media/cat/12TB/in_vivo/tim/yuki/IJ1/tif_files/IJ1pm_Apr8_stroke_30Hz.tif'
#file_name = '/media/cat/12TB/in_vivo/tim/yuki/IA3/tif_files/IA3pm_Apr8_stroke_30Hz.tif'
#file_name = '/media/cat/12TB/in_vivo/tim/yuki/IA2/tif_files/IA2pm_Apr8_stroke_30Hz.tif'
#data = tiff.imread(file_name)
#plt.imshow(np.rot90(data))
#plt.show()

#file_name = '/media/cat/12TB/in_vivo/tim/yuki/IA3/tif_files/IA3pm_Feb1_30Hz/IA3pm_Feb1_30Hz_aligned.npy'
#data = np.load(file_name)

#plt.imshow(data[1000])
#plt.show()

#np.save('/media/cat/12TB/in_vivo/tim/yuki/IA3/IA3_align_frame.npy', data[1000])

file_name = '/media/cat/12TB/in_vivo/tim/alex/vsd_june_27/VSD/31.tif'
data = tiff.imread(file_name)

file_name = '/media/cat/12TB/in_vivo/tim/alex/vsd_june_27/VSD/33.tif'
data = tiff.imread(file_name)

sum_temp = []
for k in range(len(data)):
    line = np.zeros(128, dtype=np.int32)
    line[data[k]]=1.0
    sum_temp.append(line)
        
sum_temp = np.vstack((sum_temp))
plt.imshow(sum_temp)
plt.show()
