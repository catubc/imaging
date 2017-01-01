import os
from shutil import copyfile

main_dir = '/media/cat/12TB/in_vivo/tim/dongsheng/'

cfiles = []
for root, dirs, files in os.walk(main_dir):
    for file_ in files:
        if 'time_course_data' in file_:
            src = root+'/'+ file_
            dst = src+'_barrel'
            copyfile(src, dst)
      

