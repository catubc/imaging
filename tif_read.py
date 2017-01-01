from libtiff import TIFF
import matplotlib.pyplot as plt

file_name = '/media/cat/12TB/in_vivo/tim/yuki/AQ3/tif_files/AQ3am_Apr29_Week8_30Hz.tif'

from PIL import Image

im = Image.open(file_name)

plt.imshow(im)

plt.show()
