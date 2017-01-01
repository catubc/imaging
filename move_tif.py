
#Mouse class and subclasses for processing lever pull data projec Greg Silasi and Tim Murphy lab
#Code: Cat Mitelut

import os
import glob
import numpy as np
import struct
import string, re
import scipy
import tifffile as tiff
import cPickle as pickle
import dill
import gc
from skimage.measure import block_reduce
import matplotlib.pyplot as plt
import shutil
from imreg import *
from lever_pull_utils import *




home_dir = '/media/cat/12TB/in_vivo/tim/yuki/Y1/tif_files/'
fnames = glob.glob(home_dir+"*.tif")   #use .tif extension otherwise will load .npy files

for fname in fnames: 
    if '.tif' in fname: 
        dir_name = fname.replace('.tif','/')
        os.makedirs(dir_name)

        new_name = dir_name + fname.replace(home_dir,'')
        print fname
        print new_name
        shutil.move(fname, new_name)    
    
    #quit()
