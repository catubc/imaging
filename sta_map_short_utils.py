#Utilities for sta_map_short version; basically this is a streamlined version of STA map code
import numpy as np
import multiprocessing as mp
import os
import glob
import matplotlib.pyplot as plt
from matplotlib.path import Path

class Loadptcs(object):
    """Polytrode clustered spikes file neuron record"""
    def __init__(self, sorted_file, work_dir, ptcs_flag, save_timestamps):
        
        f = open(work_dir+sorted_file+'.ptcs', "rb")
        self.sorted_file = sorted_file
        self.name = sorted_file
        self.full_path = work_dir + sorted_file
        # call the appropriate method:
        self.VER2FUNC = {1: self.readHeader, 2: self.readHeader, 3: self.readHeader}

        self.readHeader(f, ptcs_flag, save_timestamps)
        
        self.nid = []  #Make unique unit id list for loading later.
        
        self.loadData(self.nsamplebytes, f, work_dir, ptcs_flag, save_timestamps)
        
        f.close()

    def __getstate__(self):
        """Instance methods must be excluded when pickling"""
        d = self.__dict__.copy()
        try: del d['VER2FUNC']
        except KeyError: pass
        return d

    def readHeader(self, f, ptcs_flag, save_timestamps):
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


    def loadData(self, n_bytes, f, work_dir, ptcs_flag, save_timestamps):
        #call the appropriate method:
        #self.VER2FUNC = {1: self.read_ver_1, 2:self.read_ver_2, 3:self.read_ver_3}
        self.nsamplebytes = n_bytes
        self.wavedtype = {2: np.float16, 4: np.float32, 8: np.float64}[self.nsamplebytes]

        self.n_units=self.nneurons
        #print self.nneurons
        self.units=[None]*self.n_units
        self.uid = [None]*self.n_units  #Unique id for full track sorts
        self.n_sorted_spikes = [None]*self.n_units
        self.ptp=np.zeros((self.n_units), dtype=np.float32)
        self.size = []
        self.maxchan = []
        
        #print "No. of units sorted: ", self.nneurons

        #print work_dir
        for k in range(self.n_units):
            self.readUnit(f,work_dir, ptcs_flag)
            self.units[k]= self.spikes
            #print "Unit: ", k, " spikes: ", self.spikes
            if 'martin' in self.full_path:
                self.uid[k]= self.nid
            else: #All other sorts are from Nick's SS so should be the same
                self.uid[k]= self.nid-1
               
            #print "SAMPLERATE: ", self.samplerate
            if ptcs_flag: #Martin's data has wrong flag for saves
                self.units[k]=[x*self.samplerate/1E+6 for x in self.units[k]] #Converts spiketimes from usec to timesteps
            else:
                self.units[k]=[x*self.samplerate/2/1E+6 for x in self.units[k]] #Converts spiketimes from usec to timesteps

            self.n_sorted_spikes[k] = len(self.units[k])
            self.size.append(self.nspikes)
            self.maxchan.append(self.maxchanu)
            self.ptp[k]=max(self.wavedata[np.where(self.chans==self.maxchanu)[0][0]]) - \
                        min(self.wavedata[np.where(self.chans==self.maxchanu)[0][0]]) #compute PTP of template;
            #self.ptp[k]=1
            #print self.name, " unit: ", k, " ptp: ", self.ptp[k], " maxchan: ", self.maxchanu

        f.close()

        ##Save .csv spike-times files for CK
        if save_timestamps: 
            for i in range(len(self.units)):
                with open(work_dir+"timestamps_"+str(i)+".csv", "w") as f:
                    writer = csv.writer(f)
                    for j in range(len(self.units[i])):
                        writer.writerow([round(float(self.units[i][j])/float(self.samplerate)*20,6)]) #LFP events X 20

        #if (os.path.exists(work_dir+self.sorted_file+'_ptps.csv')==False):
        #    np.savetxt(work_dir+self.sorted_file+'_ptps.csv', self.ptp, delimiter=",")
        #    np.savetxt(work_dir+self.sorted_file+'_size.csv', self.size, delimiter=",")

    def readUnit(self,f, work_dir, ptcs_flag):
        self.nid = int(np.fromfile(f, dtype=np.int64, count=1)) # nid
        self.ndescrbytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # ndescrbytes
        self.descr = f.read(self.ndescrbytes).rstrip('\0 ') # descr

        if self.descr:
            try:
                self.descr = eval(self.descr) # might be a dict
            except: pass

        self.clusterscore = float(np.fromfile(f, dtype=np.float64, count=1)) # clusterscore
        self.xpos = float(np.fromfile(f, dtype=np.float64, count=1)) # xpos (um)
        self.ypos = float(np.fromfile(f, dtype=np.float64, count=1)) # ypos (um)
        self.zpos = float(np.fromfile(f, dtype=np.float64, count=1)) # zpos (um)
        self.nchans = int(np.fromfile(f, dtype=np.uint64, count=1)) # nchans
        self.chans = np.fromfile(f, dtype=np.uint64, count=self.nchans) #NB: Some errors here from older .ptcs formats
        self.maxchanu = int(np.fromfile(f, dtype=np.uint64, count=1)) # maxchanid

        self.nt = int(np.fromfile(f, dtype=np.uint64, count=1)) # nt: number of time points in template

        self.nwavedatabytes, self.wavedata = self.read_wave(f) #TEMPLATE

        self.nwavestdbytes, self.wavestd = self.read_wave(f) #STANDARD DEVIATION
        self.nspikes = int(np.fromfile(f, dtype=np.uint64, count=1)) # nspikes

        # spike timestamps (us):
        self.spikes = np.fromfile(f, dtype=np.uint64, count=self.nspikes)
        #print self.spikes
        #time.sleep(1)
        #print len(self.spikes)
        #quit()

        # convert from unsigned to signed int for calculating intervals:
        self.spikes = np.asarray(self.spikes, dtype=np.float64)

        #if 'nick' in work_dir:     #Use +1 indexes for channel #s
        if ptcs_flag: #0 = martin's sort
            #print "loading nick's data"
            #self.maxchanu-=1 #OLDER DATA
            pass
        else:
            pass
            #print "loading martin's data" #martin's data doesn't need to be converted
            #self.spikes=self.spikes #convert Martin's spike times to timesteps 
            #print self.spikes[-1]
            #quit()
            #print self.spikes
            #quit()
            
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

    def read_tsf(self,f):
        pass

def Spike_averages_parallel_prespike_3sec((args)):
    global images_temp, n_pixels, temp_window, temp_img_rate

    temp3 = args

    sum_images = np.zeros((len(temp3[0]),n_pixels,n_pixels), dtype=np.float32)
    for i in range(0,len(temp3),1):
        #baseline = np.average(images_temp[temp3[i][0:len(temp3[i])/2]], axis=0) 
        baseline = np.average(images_temp[np.arange(temp3[i][0]-temp_img_rate*3,temp3[i][0],1)], axis=0)  #Go back completely out of window for 3 sec
        
        temp_img = images_temp[temp3[i]] #Remove avg of 1st half of images
        temp_frame = (temp_img - baseline)/baseline
        sum_images += temp_frame
    
    sum_images = sum_images/max(len(temp3),1)
    
    return sum_images
    

def Compute_spike_triggered_average(unit, channel, spikes, images_raw, window, img_times,  main_dir, track_file, overwrite, img_rate, n_procs, spike_mode):
    '''Computes average frame values from t=-window .. +window (usually 180 frames for 3sec window and 30Hz sampling) '''
    
    global n_pixels, images_temp, temp_window, temp_img_rate
    
    n_pixels = len(images_raw[0])
    temp_window = window
    temp_img_rate = img_rate
    print "No. of processors: ", n_procs

    #Remove spikes outside of imaging frames
    print "Total no. spikes: ", len(spikes)
    print spikes
    print len(spikes)
    print img_times
    print len(img_times)
    
    #quit()
    
    temp0 = np.where(np.logical_and(spikes>=img_times[0]+2*window, spikes<=img_times[-1]-window))[0]    #Exclude spikes too close to beginning or end of recordings.
    spikes = spikes[temp0]

    temp3 = []
    if spike_mode =='all': #Use all spikes
        for spike in spikes:        #temp3 contains all the frame indexes from img_times for each spike in raster; e.g. 180 frames for each spike automatically aligned
            temp3.append(np.where(np.logical_and(img_times>=spike-window, img_times<=spike+window))[0][0:2*int(window*img_rate)]) #Fixed this value or could be off +/-1 frame
    elif spike_mode =='lockout': 
        for s in range(1, len(spikes)-1,1):        #temp3 contains all the frame indexes from img_times for each spike in raster; e.g. 180 frames for each spike automatically aligned
            if ((spikes[s]-window)>=spikes[s-1]) and ((spikes[s]+window)<=spikes[s+1]):
                array_temp = np.where(np.logical_and(img_times>=spikes[s]-window, img_times<=spikes[s]+window))[0][0:2*int(window*img_rate)]
                temp3.append(array_temp) #Fixed this value or could be off +/-1 frame
    elif spike_mode =='burst': 
        for s in range(1, len(spikes),1):        #temp3 contains all the frame indexes from img_times for each spike in raster; e.g. 180 frames for each spike automatically aligned
            if ((spikes[s]-window)>=spikes[s-1]):
                array_temp = np.where(np.logical_and(img_times>=spikes[s]-window, img_times<=spikes[s]+window))[0][0:2*int(window*img_rate)]
                temp3.append(array_temp) #Fixed this value or could be off +/-1 frame

    print "No. spks in img period: ", len(spikes), " rate: ", float(len(spikes))/(img_times[-1]-img_times[0]), " Hz, no. spks pass criteria: ", len(temp3)
        
    if len(temp3)==0: return 

    #Check to see if images already loaded and saved as .npy
    #stm_file_name = main_dir + 'stm_files/' + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_'+plot_string+'_'+str(window)+'sec_window_'+str(len(spikes)).zfill(5)+"_spikes"
    stm_file_name = main_dir + 'stm_files/img_avg_' + track_file+'_unit'+str(unit).zfill(3)+'_ch'+str(channel).zfill(3)+'_'+spike_mode+'_'+str(window)+'sec_window_'
    
    if (overwrite) or (len(glob.glob(stm_file_name+"*"))==0):
    #if overwrite :
        #images_temp = np.array(images_raw.copy(), dtype=np.float32)
        #images_temp = images_raw.astype(np.float32, copy=False)
        images_temp = images_raw
        
        #Remove baseline - global activity regression
        #if False:
        #    print "Removing baseline over all data"
        #    baseline = np.mean(images_temp, axis=0)
        #    #baseline = np.mean(images_temp, axis=0, dtype=np.float32)
        #    images_temp = (images_temp - baseline)/baseline
        #else: print "Not removing baseline"

        #Compute all frames based on image index;
        print "... stm processing in parallel for window: ", window, " secs ..."
        images_triggered_temp=[]
        temp4 = []  #Contains only first spikes from bursts in cell

        #*********** USE PARALLEL CODE ************
        if False: 
            if len(temp3) < 30:
                pool = mp.Pool(1)
                temp4.append(temp3)            
            else:
                pool = mp.Pool(n_procs)
                chunks = int(len(temp3)/n_procs) #Break up the temp3 array into n_procs that are "chunk" long each
                for i in range(n_procs):
                    temp4.append(temp3[i*chunks:(i+1)*chunks])
                    
                #DISCARD RESIDUE: May wish to recapture this eventually.
                #if n_procs*chunks<len(spikes):
                #    temp_sum = images_temp[temp3[n_procs*chunks]]
                #    for i in range(n_procs*chunks+1,len(spikes),1):
                #        temp_sum+=images_temp[temp3[i]]
                #    images_triggered_temp.append(temp_sum)
                
            #print "Removing average of all pre spike frames - (time: -", window, "sec .. 0sec)"
            images_triggered_temp.extend(pool.map(Spike_averages_parallel_prespike_3sec, temp4))
            
            pool.close()
            print "... done "

            #Sum over all spikes
            print "Summing Number of chunks: ", len(images_triggered_temp)

            temp_images = np.zeros((int(window*img_rate)*2, n_pixels, n_pixels), dtype=np.float16)
            for i in range(len(images_triggered_temp)):
                temp_images += images_triggered_temp[i]
            
            #DIVIDE BY NUMBER OF CHUNKS; Note used to be divided by number of spikes; also residue is being thrown out...
            images_processed = temp_images/float(len(images_triggered_temp))

        #************ USE SINGLE CORE TO COMPUTE AVERAGES **************
        else:
            temp3 = np.array(temp3)
            images_processed = np.zeros((len(temp3[0]),n_pixels,n_pixels), dtype=np.float32)
            for i in range(0,len(temp3),1):
                print "...averaging spike: ", i, " / ", len(temp3)
                #baseline = np.average(images_temp[temp3[i][0:len(temp3[i])/2]], axis=0) 
                baseline = np.average(images_temp[np.arange(temp3[i][0]-int(temp_img_rate)*3,temp3[i][0],1)], axis=0)  #Go back 3 sec prior to temp3[0] value...
                
                temp_img = images_temp[temp3[i]] #Remove avg of 1st half of images
                temp_frame = (temp_img - baseline)/baseline
                images_processed += temp_frame
            
            images_processed = images_processed/max(len(temp3),1)
                
                
        #npy_file_name = file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_'+plot_string+'_'+str(window)+'sec_window_'+str(len(spikes)).zfill(5)+"_spikes"
        stm_file_name= stm_file_name+str(len(temp3)).zfill(5)+"_spikes"

        np.save(stm_file_name, images_processed)

    else: 
        print "Skipping processing of images ... loading from file"
        images_processed = np.load(glob.glob(stm_file_name+"*")[0])

    print "Shape mean: ", images_processed.shape
    
def Define_generic_mask(images_processed, main_dir):

    global coords, images_temp, ax, fig, cid
    
    images_temp = images_processed
    
    fig, ax = plt.subplots()

    if (os.path.exists(main_dir + 'genericmask.txt')==False):
        coords=[]

        ax.imshow(images_processed[100])#, vmin=0.0, vmax=0.02)
        ax.set_title("Compute generic (outside the brain) mask")
        #figManager = plt.get_current_fig_manager()
        #figManager.window.showMaximized()
        cid = fig.canvas.mpl_connect('button_press_event', on_click)
        plt.show()

        #******* MASK AND DISPLAY AREAS OUTSIDE GENERAL MASK 
        #Search points outside and black them out:
        all_points = []
        for i in range(len(images_processed[0][0])):
            for j in range(len(images_processed[0][0])):
                all_points.append([i,j])

        all_points = np.array(all_points)
        vertixes = np.array(coords) 
        vertixes_path = Path(vertixes)
        
        mask = vertixes_path.contains_points(all_points)
        counter=0
        coords_save=[]
        for i in range(len(images_processed[0][0])):
            for j in range(len(images_processed[0][0])):
                if mask[counter] == False:
                    images_processed[100][i][j]=0
                    coords_save.append([i,j])
                counter+=1

        fig, ax = plt.subplots()
        ax.imshow(images_processed[100])
        plt.show()
       
        genericmask_file = main_dir + 'genericmask.txt'
        np.savetxt(genericmask_file, coords_save)

        print "Finished Making General Mask"

    #else:
    #    print "Loading saved general mask"
        
        
    #if (os.path.exists(main_dir + 'bregmamask.txt')==False):
        #bregma_coords = []
        #print "Making Bregma mask"
        #ax.imshow(images_processed[100])#, vmin=0.0, vmax=0.02)
        #ax.set_title("Compute bregma mask")
        ##figManager = plt.get_current_fig_manager()
        ##figManager.window.showMaximized()
        #cid = fig.canvas.mpl_connect('button_press_event', remove_bregma)
        #plt.show()

       
        #bregmamask_file = main_dir + 'bregmamask.txt'
        #np.savetxt(bregmamask_file, bregma_coords)

        #print "Finished Bregma Mask"

    #else:
    #    print "Loading saved bregma mask"
        
    return generic_coords

def on_click(event):
    
    global coords, images_temp, ax, fig, cid
    
    n_pix = len(images_temp[0])
    
    if event.inaxes is not None:
        coords.append((event.ydata, event.xdata))
        for j in range(len(coords)):
            for k in range(3):
                for l in range(3):
                    images_temp[100][min(n_pix,int(coords[j][0])-1+k)][min(n_pix,int(coords[j][1])-1+l)]=0

        ax.imshow(images_temp[100])
        #plt.show()
        fig.canvas.draw()
                    #figManager = plt.get_current_fig_manager()
                    #figManager.window.showMaximized()
    else:
        print 'Exiting'
        plt.close()
        fig.canvas.mpl_disconnect(cid)
        
