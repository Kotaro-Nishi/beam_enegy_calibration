import cv2
import sys,os
from natsort import natsorted

LengthOfData = 825

arg = sys.argv

tiflist = []

k0=0

DataPath = arg[1]

filelist = os.listdir(arg[1])
for i in range(len(filelist)):
    if filelist[i][-4:] == "tiff" or filelist[i][-3:] == "tif":
        tiflist.append(filelist[i])
filelist = natsorted(tiflist)



for ii in range(k0,k0+LengthOfData):
    for i in range(4):
        im = cv2.imread(DataPath + "/" + filelist[i+ii*4],cv2.IMREAD_UNCHANGED)
        #try:
        #    l = len(im)
        #    print(l)
        #except:
        #    print(i+ii*4)
        print(type(im))