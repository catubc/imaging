import numpy as np


np.set_printoptions(precision=4)
np.set_printoptions(suppress=True)

data = np.sort(np.random.random(1000))*581

print data

file_out = '/media/cat/8TB/in_vivo/tim/dongsheng/2015-7-22/2015-7-22-5/unit_27_channel_16_ptp_050.csv'

np.savetxt(file_out, data)
