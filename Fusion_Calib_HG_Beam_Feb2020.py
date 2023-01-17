################################################################################################################
#
# Plot horizontal slices
#
#################################################################################################################

import matplotlib.pyplot as plt
import numpy as np
import cv2
import os
import csv
from scipy import misc
import glob

#from os import listdir
from scipy.optimize import curve_fit, minimize_scalar
from scipy import asarray as ar,exp

from matplotlib import rcParams
from natsort import natsorted, ns

def gaus(x,a,x0,sigma):
	return a*exp(-(x-x0)**2/(2*sigma**2))

def TwoGaus(x,a1,a2,x01,x02,sigma1,sigma2,b):
	return a1*exp(-(x-x01)**2/(2*sigma1**2))+a2*exp(-(x-x02)**2/(2*sigma2**2))+b

#rcParams['text.usetex'] = True
#rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

filelist = []
tiflist = []
mypixy = []

Commonpath = str(os.path.dirname(os.path.abspath(__file__)))

print(os.listdir(Commonpath))
Picturepath = input("Drag and drop file from folder: %s\ + ? "%(Commonpath))

img = cv2.imread(Picturepath,cv2.IMREAD_UNCHANGED)   # Load an color image in grayscale
print(img[0])
img = np.float32(img)
extends=[]
try:
    with open("extends.csv","r") as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        for row in spamreader:
            extends.append(row)

    extends = [int(extends[0][0]),int(extends[0][1])]

    print(extends[0])
    print(extends[1])
except IOError:
    print("no vertical size defined such file must be created first!")
#seems to work better if sigma presumed initially smaller than expected
#working p_initial_Four = [11000,1500,17000,27000, 130,140,333,444, 5,1,5,5]

while True:
    option = input("Obtain dispersion and 404 location: d\n Obtain vertical size light: v\n ")
    if option == "d":
        scalex =  1.
        scaley = .02
        img = cv2.resize(img,None,fx=scalex, fy=scaley, interpolation = cv2.INTER_AREA)
        DispersionList = []
        list404nm = []
        for rowpix in range(int(int(extends[0])*scaley),int(int(extends[1])*scaley)):

            x = range(len(img[rowpix]))
            y = img[rowpix]

            #p_initial_Two = [1.01,1.01, 1315., 1964., 1,.5, 100.] 20200319
            p_initial_Two = [1.01,1.01, 1238., 1892., 1,.5, 100.]
            popt,pcov = curve_fit(TwoGaus,x,y, p0 = p_initial_Two)
            print(popt)
            print(np.sqrt(pcov[0][0]), ",",np.sqrt(pcov[1][1]), ",",np.sqrt(pcov[2][2]), ",",np.sqrt(pcov[3][3]))

            DispersionList.append((407.7837-404.6565)/(popt[3]-popt[2]))
            list404nm.append(popt[2])
            print(DispersionList[-1])
            
            print("next",rowpix)
        dispersion = np.mean(DispersionList)
        dispersionsigma = np.std(DispersionList,ddof=1)
        pixelof404nm = np.mean(list404nm)
        list404nmsigma = np.std(list404nm,ddof=1)
        
        print("mean dispersion: " + str(dispersion))
        print("sigma dispersion: " + str(dispersionsigma)) #ddof=1 because n-1

        print("mean 404: " + str(np.mean(list404nm)))
        print("sigma 404: " + str(list404nmsigma)) #ddof=1 because n-1

        xfit = np.linspace(x[0], x[len(x)-1], num=len(x)*10)
        yfit = TwoGaus(xfit,*popt)

        fig5 = plt.figure()

        ax = fig5.add_subplot(111)
        ax.xaxis.grid(True)
        ax.yaxis.grid(True)

        ax.plot(x,y,'r',lw=1.3)
        ax.plot(xfit,yfit,'g',lw=1.3)
        plt.show()
        answer1 = input("Write to file y/n?\n")
        if answer1 == "y":
            with open("specdata.csv", 'w', newline='') as write_obj:
                # Create a writer object from csv module
                csv_writer = csv.writer(write_obj)
                # Add contents of list as last row in the csv file
                csv_writer.writerows([[dispersion],[dispersionsigma],[pixelof404nm],[list404nmsigma]])

    if option == "v":
        scalex =  .1
        scaley = 1.
        img_v = cv2.resize(img,None,fx=scalex, fy=scaley, interpolation = cv2.INTER_AREA)
        stdofy=[]
        pixelsinthecone=[]
        for rowpix in range(len(img_v)):
            y = img_v[rowpix]
            stdofy.append(np.std(y,ddof=1))
        threshold = max(stdofy) - (max(stdofy)-min(stdofy))*.3
        for i in range(len(stdofy)):
            if stdofy[i] > threshold:
                pixelsinthecone.append(i)
        upperborder = min(pixelsinthecone)
        lowerborder = max(pixelsinthecone)
        print("upper border pixel: " + str(upperborder))
        print("lower border pixel: " + str(lowerborder))
        fig5 = plt.figure()

        ax = fig5.add_subplot(111)
        ax.xaxis.grid(True)
        ax.yaxis.grid(True)

        ax.plot(range(len(stdofy)),stdofy,'r',lw=1.3)

        plt.show()
        #print(type(upperborder))
        answer2 = input("Write to file y/n?\n")
        if answer2 == "y":
            extends = [upperborder,lowerborder]
            with open("extends.csv", 'w', newline='') as write_obj:
                # Create a writer object from csv module
                csv_writer = csv.writer(write_obj)
                # Add contents of list as last row in the csv file
                csv_writer.writerow([upperborder,lowerborder])