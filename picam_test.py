import numpy as np
import matplotlib.pyplot as plt

#img_filename = '/media/pi/ssd/lolcat_128.raw'

img_filename = '/media/cat/8TB/picam/lolcat_128_ram.raw'
f = open(img_filename, 'rb')

data = f.read()
print data[:10].decode('utf-8')
quit()


filename = '/media/pi/Cage_4A/lolcat_128_ram.raw_time.txt'
filename = '/media/cat/8TB/picam/lolcat_128_ram_60Hz.raw_time.txt'
#filename = '/media/cat/8TB/picam/lolcat_128_ram.raw_time.txt'

bin_width = 0.00001

data = np.loadtxt(filename)
print data
data = data*1E-6
plt.plot(data)
plt.show()

isi = data[1:]-data[:-1]

if False: 
    print "...# frames: ", len(data)
    indexes = np.where(isi>0.035)[0]
    print "... time of skipped frames: ", data[indexes]
    print "... isis skipped (secs): ", isi[isi>0.035]
    print "... total skpped time (secs): ", np.sum(isi[isi>0.035])
    y = np.histogram(isi, bins = np.arange(0.030,0.040,bin_width))
else: #60Hz
    print "...# frames: ", len(data)
    indexes = np.where(isi>0.017)[0]
    print "... time of skipped frames: ", data[indexes]
    print "... isis skipped (secs): ", isi[isi>0.017]
    print "... total skpped time (secs): ", np.sum(isi[isi>0.017])
    y = np.histogram(isi, bins = np.arange(0.012,0.020,bin_width))

plt.bar(y[1][:-1], y[0], bin_width, color='blue')
plt.ylim(0,10)
plt.show()
