import matplotlib.pyplot as plt

from numpy import transpose as nptp
import numpy as np

import cv2
import os , sys
import re

import csv

from natsort import natsorted
from scipy.optimize import curve_fit
from lmfit import minimize, Parameters, report_fit, Minimizer

import random

def sinefit(x,a,gamma,phi,e,b):
	return a * np.sin(2*np.pi*x/(2.*WLofPixel*gamma**2) + phi)*(1+e*x) + b

NumOfStripes = 1
NumOfSeries = 825
LengthOfInterval = 825

dispersion = 0.004808759252942803e-9 #Former 0.0048085643228767475e-9
pixelof404nm = 1315.3024745554976 #Former: 1243.8551272671261
CalibWl = 404.6565e-9

#(407.7837-404.6565)

chosenpx = 1300

global listOfboots
listOfboots = []

LambdaRfit = 412e-9

arg = sys.argv

tiflist = []

filelist = os.listdir(arg[1])
for i in range(len(filelist)):
    if filelist[i][-4:] == "tiff" or filelist[i][-3:] == "tif":
        tiflist.append(filelist[i])
filelist = natsorted(tiflist)

Commonpath = str(os.path.dirname(os.path.abspath(__file__)))

print(Commonpath)

DataPath = arg[1]

print(DataPath)

slashes = [i for i,ltr in enumerate(DataPath) if ltr == "\\"]
DataPathBase = DataPath[0:slashes[-2]]

print(DataPathBase)

Backgroundimage = cv2.imread(DataPathBase + "/BG/raw_shortExp/File1.tif",cv2.IMREAD_UNCHANGED)
ListOfAllMatrices = []
ListOfAllSigmaMatrices = []
print(DataPathBase + "/BG/raw_shortExp/File1.tif")

WidthOfStripe = 5
print(WidthOfStripe)

MatrixOfSlit = []

maybeBadGamma = 0

badk = []

DoNotSave = 1

FirstWL = (0 - pixelof404nm)*dispersion + CalibWl
WLofPixel = (chosenpx - pixelof404nm)*dispersion + CalibWl
LastWL = (2304 - pixelof404nm)*dispersion + CalibWl

print("First wavelength: ", FirstWL)
print("Wavelength at which the analysis is done: ", WLofPixel)
print("Last wavelength: ", LastWL)

k0 = 000
k = k0
SeriesListPictureMinusBGWidthOfStripe4Mean=[]

for ii in range(k0,k0+LengthOfInterval):
    ListPictureMinusBGWidthOfStripe = []
    print(filelist[ii*4])
    for i in range(4):
        noti = i
        PictureMinusBG = cv2.imread(DataPath + "/" + filelist[i+ii*4],cv2.IMREAD_UNCHANGED).astype(int)-Backgroundimage

        PictureMinusBG = np.float32(PictureMinusBG)
        
        PictureMinusBGWidthOfStripe = nptp(nptp(PictureMinusBG)[(chosenpx-int(WidthOfStripe/2)):chosenpx+int(WidthOfStripe/2)+1])
        
        ListPictureMinusBGWidthOfStripe.append(PictureMinusBGWidthOfStripe)

    ListPictureMinusBGWidthOfStripe4Mean = sum(ListPictureMinusBGWidthOfStripe)/4.
    
    SeriesListPictureMinusBGWidthOfStripe4Mean.append(ListPictureMinusBGWidthOfStripe4Mean)

print(NumOfSeries)
checkstart = 1

residsqrlst=[]

for Saves in range(0,1):

    print("Saves " + str(Saves))

    dataBooted = []
    
    residLst = []

    interpolatedpicture = []

    for ii in range(LengthOfInterval):
        stripe = SeriesListPictureMinusBGWidthOfStripe4Mean[ii]

        stripe = cv2.resize(stripe,(1,len(ListPictureMinusBGWidthOfStripe4Mean)),interpolation=cv2.INTER_AREA) #1 zu 1000

        stripe = nptp(stripe)[0]

        stripejjj=[]
        for jjj in range(len(stripe)):
            stripejjj.append(stripe[jjj])
            
        stripe = nptp(np.array([stripejjj]))
        
        interpolatedpicture.append(nptp(stripe)[0])
        
    data = np.array(interpolatedpicture)
    """
    for Boot in range(NumOfSeries):
        dataBooted.append(data[listOfboots[Boot]])
    dataBooted = np.array(dataBooted)

    dataLst = []
    for Strpidx in range(NumOfSeries):
        dataLst.append(dataBooted[Strpidx][500])
    """
    dataLst = nptp(data)[510]
    
    xlst = np.array(range(NumOfSeries))/1000.

    fntsze = 13
    fig = plt.figure("Sum of squared resuduals vs. Position on camera screen")
    ax = fig.add_subplot(111)
    ax.set_title("Sum of squared residuals vs. Position on camera screen")
    ax.set_xlabel(r"Position on camera screen [$\mu$m]",fontsize=fntsze)
    ax.set_ylabel(r"Sum of squared residuals [a.u.]",fontsize=fntsze)
    ax.plot(xlst,dataLst)
    plt.show()
    
    Guessvalues = [100.4,352.4,.2,0.,100.]#,.35e-4]

    popt, pcov = curve_fit(sinefit,xlst,dataLst,Guessvalues)
    
    print(popt)
    
    fntsze = 13
    fig = plt.figure("Sum of squared resuduals vs. Position on camera screen")
    ax = fig.add_subplot(111)
    ax.set_title("Sum of squared residuals vs. Position on camera screen")
    ax.set_xlabel(r"Position on camera screen [$\mu$m]",fontsize=fntsze)
    ax.set_ylabel(r"Sum of squared residuals [a.u.]",fontsize=fntsze)
    ax.plot(xlst,dataLst)
    ax.plot(xlst,sinefit(xlst,*popt))
    ax.plot(xlst,sinefit(xlst,Guessvalues[0],Guessvalues[1],Guessvalues[2],Guessvalues[3],Guessvalues[4]))
    plt.show()
    
    #residsqrlst.append(np.sum(residLst**2))

    #UnDi.SingleAmplitudeDone = 0

print(residsqrlst)
fntsze = 13
fig = plt.figure("Sum of squared resuduals vs. Position on camera screen")
ax = fig.add_subplot(111)
ax.set_title("Sum of squared residuals vs. Position on camera screen")
ax.set_xlabel(r"Position on camera screen [$\mu$m]",fontsize=fntsze)
ax.set_ylabel(r"Sum of squared residuals [a.u.]",fontsize=fntsze)
ax.plot(np.array(range(len(residsqrlst)))*10.+190.,residsqrlst)
plt.show()