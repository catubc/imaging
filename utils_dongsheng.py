import struct, array, csv
from tsf_ptcs_classes import *


def Low_pass_tsf(tsf):

    #probe = [16, 1, 15, 2, 12, 5, 13, 4, 10, 7, 9, 8, 11, 6, 14, 3]
        
    file_name = tsf.tsf_name[:-8] + '_lp.tsf'
    f1 = open(file_name, 'wb')

    f1.write(tsf.header)
    f1.write(struct.pack('i', tsf.iformat))
    f1.write(struct.pack('i', tsf.SampleFrequency))
    f1.write(struct.pack('i', tsf.n_electrodes))
    f1.write(struct.pack('i', tsf.n_vd_samples))
    f1.write(struct.pack('f', tsf.vscale_HP))
    print tsf.Siteloc
    for i in range (tsf.n_electrodes):
        f1.write(struct.pack('h', tsf.Siteloc[i*2]))
        f1.write(struct.pack('h', tsf.Siteloc[i*2+1]))
        f1.write(struct.pack('i', i+1)) #Need to add extra value for Fortran arrays

    ecp_lp=tsf.ec_traces.copy()
    ecp_temp = tsf.ec_traces.copy()

    ecp_lp = ecp_temp - wavelet(ecp_lp, wname="db4", maxlevel=6)
    
    ecp_temp = []
    for i in range(len(ecp_lp)):
        ecp_temp.append(ecp_lp[i][:tsf.n_vd_samples])
    ecp_temp = np.int16(ecp_temp)
    
    print "Writing data"
    for i in range(len(ecp_lp)):
        ecp_temp[i].tofile(f1)
        
        #ecp_temp[probe[i]-1][::20].tofile(f1) #USE THIS FORMULA TO SUBSAMPLE IF REQUIRED

    f1.write(struct.pack('i', 0)) #Write # of fake spikes
    f1.close()
    print "DONE"
    
def Classify_traces(location, hemisphere, maxormin, state_match, tc_data, area_maxormin, area_hemisphere, depth_left, state_left, 
                    channel_left, raster_left, template_left, area_index_left, mouse_index, file_index):
    
    '''Splits all data into 4 classes: depth x state. Outputs array of selected class back to main routine
    '''

    cortex_traces_awake = []
    cortex_traces_anesthetized = []
    subcortical_traces_awake = []
    subcortical_traces_anesthetized = []

    cortex_awake_channels=[]
    cortex_anesthetized_channels=[]
    subcortical_awake_channels=[]
    subcortical_anesthetized_channels=[]

    cortex_awake_rasters=[]
    cortex_anesthetized_rasters=[]
    subcortical_awake_rasters=[]
    subcortical_anesthetized_rasters=[]
    
    cortex_awake_templates=[]
    cortex_anesthetized_templates=[]
    subcortical_awake_templates=[]
    subcortical_anesthetized_templates=[]

    cortex_awake_mouse=[]
    cortex_anesthetized_mouse=[]
    subcortical_awake_mouse=[]
    subcortical_anesthetized_mouse=[]

    cortex_awake_file_indexes_out = []
    cortex_anesthetized_file_indexes_out = []
    subcortical_awake_file_indexes_out = []
    subcortical_anesthetized_file_indexes_out = []
    
    
    cortical_states = ['awake','anesthetized','awake','anesthetized']
    cortical_depths = ['cortex','cortex','anesthetized','anesthetized']
    array_max = []
    depth=[]
    state=[]
    area_index = location      
    #area_names = [ 'hindlimb', 'forelimb', 'barrel', 'motor', 'visual', 'retrosplenial', 'acc', 'allcortex']
     
    #CLASSIFY CELLS BY DEPTH AND ANESTHETIC
    print "Looking for: ", location, hemisphere, cortical_states[state_match]
    print "Cells in...: ", len(tc_data)
    
    from collections import Counter
    print Counter(state_left)

    #print "Awake cells...", state_left.count['awake']
    
    #counter=0
    for i in range(len(tc_data)):
        if max(tc_data[i])==min(tc_data[i]): continue #SKIP courses set to = 0 in main routine (i.e. missed thresholds)
        print "good cell"
        #if (area_index_left[i] != location) or (area_hemisphere[i] != hemisphere) or (area_maxormin[i] != maxormin): continue
        if (area_index_left[i] != location) or (area_hemisphere[i] != hemisphere): continue
        
        if (area_maxormin[i]!='min'): continue #Start on "min" and load next two traces for 'min' and 'max'

        #Load the trace with largest ABS value from -1sec .. +3sec
        if max(abs(tc_data[i][200:]))>max(abs(tc_data[i+1][200:])):
            array_max = tc_data[i]  #area_location, hemisphere, states,
        else:
            array_max = tc_data[i+1]  #area_location, hemisphere, states,
            
        depth = depth_left[i]
        state = state_left[i]
        mouse = mouse_index[i]
            
        if depth=='cortex':
            if state=='awake':
                cortex_traces_awake.append(array_max)
                cortex_awake_channels.append(channel_left[i])
                cortex_awake_rasters.append(raster_left[i])
                cortex_awake_templates.append(template_left[i])
                cortex_awake_mouse.append(mouse)
                cortex_awake_file_indexes_out.append(file_index[i])
                #print file_index[i]
            else:
                cortex_traces_anesthetized.append(array_max)
                cortex_anesthetized_channels.append(channel_left[i])
                cortex_anesthetized_rasters.append(raster_left[i])
                cortex_anesthetized_templates.append(template_left[i])
                cortex_anesthetized_mouse.append(mouse)
                cortex_anesthetized_file_indexes_out.append(file_index[i])
        else:   
            if state=='awake':
                subcortical_traces_awake.append(array_max)
                subcortical_awake_channels.append(channel_left[i])
                subcortical_awake_rasters.append(raster_left[i])
                subcortical_awake_templates.append(template_left[i])
                subcortical_awake_mouse.append(mouse)
                subcortical_awake_file_indexes_out.append(file_index[i])

            else:
                subcortical_traces_anesthetized.append(array_max)
                subcortical_anesthetized_channels.append(channel_left[i])
                subcortical_anesthetized_rasters.append(raster_left[i])
                subcortical_anesthetized_templates.append(template_left[i])
                subcortical_anesthetized_mouse.append(mouse)
                subcortical_anesthetized_file_indexes_out.append(file_index[i])
    
    cluster_traces_array = []
    cluster_traces_array.append(cortex_traces_awake)
    cluster_traces_array.append(cortex_traces_anesthetized)
    cluster_traces_array.append(subcortical_traces_awake)
    cluster_traces_array.append(subcortical_traces_anesthetized)
    
    channels=[]
    channels.append(cortex_awake_channels)
    channels.append(cortex_anesthetized_channels)
    channels.append(subcortical_awake_channels)
    channels.append(subcortical_anesthetized_channels)

    rasters=[]
    rasters.append(cortex_awake_rasters)
    rasters.append(cortex_anesthetized_rasters)
    rasters.append(subcortical_awake_rasters)
    rasters.append(subcortical_anesthetized_rasters)

    templates=[]
    templates.append(cortex_awake_templates)
    templates.append(cortex_anesthetized_templates)
    templates.append(subcortical_awake_templates)
    templates.append(subcortical_anesthetized_templates)
    
    mouse = []
    mouse.append(cortex_awake_mouse)
    mouse.append(cortex_anesthetized_mouse)
    mouse.append(subcortical_awake_mouse)
    mouse.append(subcortical_anesthetized_mouse)
    
    files = []
    files.append(cortex_awake_file_indexes_out)
    files.append(cortex_anesthetized_file_indexes_out)
    files.append(subcortical_awake_file_indexes_out)
    files.append(subcortical_anesthetized_file_indexes_out)

    #print "... cells out: ", len(cluster_traces_array[2])

    #window = 3
    #n_interp_points = 600
    #xx = np.linspace(-float(window),float(window),n_interp_points)

    #ave_trace = np.zeros(len(xx),dtype=np.float32)
    #for k in range(len(cluster_traces_array[state_match])):
        #plt.plot(xx,np.array(cluster_traces_array[state_match]), color='magenta', linewidth=10, alpha=1)
        #ave_trace+=cortex_traces_awake[k]
        #plt.plot([0,0],[-3,3], color='black', linewidth=2)
        #plt.plot([-3,3],[0,0], color='black', linewidth=2)
        #plt.ylim(-1,3)
        #plt.xlim(-3,3)
        #plt.tick_params(axis='both', which='major', labelsize=30)
        #plt.yticks([-1,0,1,2,3])
    #plt.show()
    ##plt.plot(xx, ave_trace/(k-1), color='black', linewidth=10,alpha=1)
    ##plt.show()
    
    print "... cells out: ", len(cluster_traces_array[state_match])
    print ".."
    
    print "NO. filenames: ", len(files[state_match])
    
    return cluster_traces_array[state_match], channels[state_match], rasters[state_match], templates[state_match], mouse[state_match],files[state_match]

def draw_pie(ax,ratios=[0.4,0.3,0.3], X=0, Y=0, size=1000, colors=colors):
    
    #colors=['black','brown','violet','dodgerblue','mediumvioletred','indianred','lightseagreen','lightsalmon','pink','darkolivegreen']

    N = len(ratios)
 
    xy = []
 
    start = 0.
    for ratio in ratios:
        x = [0] + np.cos(np.linspace(2*math.pi*start,2*math.pi*(start+ratio), 30)).tolist()
        y = [0] + np.sin(np.linspace(2*math.pi*start,2*math.pi*(start+ratio), 30)).tolist()
        xy1 = zip(x,y)
        xy.append(xy1)
        start += ratio
 
    for i, xyi in enumerate(xy):
        ax.scatter([X],[Y] , marker=(xyi,0), s=size, facecolor=colors[i])
