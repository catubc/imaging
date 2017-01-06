from matplotlib import pyplot as plt
import numpy as np
import scipy
import scipy.misc
import struct

filename = '/media/cat/8TB/picam/old/lolcat_shell.bin'
fin = open(filename, "rb")
a = np.fromfile(fin, dtype=np.uint32)*1E-6
a = a[a>0]  #Remove zeros at end of array; 

plt.plot(a)
plt.show()

##filename = '/media/cat/ssd/lolcat_1.raw_time.txt'
#filename = '/media/cat/8TB/picam/lolcat_2.raw_time.txt'
#a = np.loadtxt(filename)
#a = a-a[0]

#compute isis
isi = a[1:]-a[:-1]

print "... # frames: ", len(a)
print "... isi skipped (secs): ", isi[isi>0.04]
print "... total skipped time (secs): ", np.sum(isi[isi>0.04])

#plt.plot(a)
#plt.show()


#quit()

bin_width = 0.001   # histogram bin width in usec
y = np.histogram(isi, bins = np.arange(0,1,bin_width))
plt.bar(y[1][:-1], y[0], bin_width, color='blue')
plt.show()
