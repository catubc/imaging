from matplotlib import pyplot as plt

import numpy as np

file_name = '/home/cat/Downloads/cat.txt'

data = np.loadtxt(file_name)
data = data[:,1]
print "...Sample inter-frame intervals (sec): ", 1./data[:5]

time_dropped = 0
total_time = 0
for k in range(len(data)):
    total_time+=1./data[k]
    
    if data[k]<28: 
        print "...error frame: ", k, " length of frame: ", 1./data[k], "sec."
        time_dropped+=1./data[k]


print "Time dropped: ", time_dropped, "sec."
print "Total time: ", total_time, "sec."
print "Expected time based on 28.81Hz x 5000 frames: ", 0.034705*5000, "sec."
