import pylab
import matplotlib.pyplot as plt
import imageio
filename = '/media/cat/12TB/in_vivo/tim/yuki/AQ2/AQ2am_Jan18_30Hz_reward_code_02.mp4'
filename2 = '/media/cat/12TB/in_vivo/tim/yuki/AQ2/AQ2am_Jan18_30Hz_reward_code_04.mp4'
filename3 = '/media/cat/12TB/in_vivo/tim/yuki/AQ2/AQ2am_Jan18_30Hz_reward_code_07.mp4'

vid = imageio.get_reader(filename,  'ffmpeg')
all_data = vid.get_data()
print all_data.shape
quit()

plt.imshow(vid.get_data(0))
plt.show()

for num in nums:
    image = vid.get_data(num)
    fig = pylab.figure()
    fig.suptitle('image #{}'.format(num), fontsize=20)
    pylab.imshow(image)
pylab.show()
