import matplotlib.pyplot as plt
import numpy as np
from numpy import transpose as nptp
import csv
import os, sys
from natsort import natsorted, ns
from scipy.optimize import curve_fit
import pandas as pd

arg = sys.argv

from scipy.spatial import distance

import random
import statistics

dp=1

chosenfolders = [
                "SharedP_20210401_825_825_seed10000_new_zpfit_boot000425",
                "SharedP_20210401_825_825_seed10000_new_zpfit_boot000825",
                "SharedP_20210401_825_825_seed10000_new_zpfit_boot400825",
                "SharedP_20210401_825_825_seed10000_new_zpfit_delta0fit",
                "SharedP_20210401_825_825_seed10000_new_zpfit_MultiWL",
                "SharedP_20210401_825_825_seed10000_new_zpfit_Yscreen"
                ]

targetfiles = ["Amplfit.csv",
               "Ytran.csv",
               "gammafit.csv",
               "delta0fit.csv",
               "cfit.csv",
               "Amplfitstderr.csv",
               "Ytranstderr.csv",
               "gammafitstderr.csv",
               "delta0fitstderr.csv",
               "cfitstderr.csv",
               "amplYtranCorrel.csv",
               "amplgammafitCorrel.csv",
               "ampldelta0fitCorrel.csv",
               "amplcfitCorrel.csv",
               "YtrangammafitCorrel.csv",
               "Ytrandelta0fitCorrel.csv",
               "YtrancfitCorrel.csv",
               "gammadelta0fitCorrel.csv",
               "gammacfitCorrel.csv",
               "delta0fitcfitCorrel.csv",
               "zpfit.csv",
               "Residualmean.csv"]

def h2(xdata,ydata):
    x_min = np.min(xdata)
    x_max = np.max(xdata)
      
    y_min = np.min(ydata)
    y_max = np.max(ydata)
      
    x_bins = np.linspace(x_min, x_max, 60)
    y_bins = np.linspace(y_min, y_max, 60)
      
    fig, ax = plt.subplots(figsize =(10, 10))
    # Creating plot
    plt.hist2d(xdata, ydata, bins =[x_bins, y_bins])
    plt.title("Changing the bin scale")
      
    ax.set_xlabel('X-axis') 
    ax.set_ylabel('X-axis') 
      
    # show plot
    plt.tight_layout() 
    if dp:
        plt.show()
    else:
        print("dp is set to 0 use plt.show()")

def gaussfitfunc(x, amp, ev, sigma):
	return amp/np.sqrt(2*np.pi*sigma**2) * np.exp(-(x-ev)**2/(2*sigma**2))

def gaussfitfunc2(x, amp1, ev1, sigma1, amp2, ev2, sigma2):
	return amp1/np.sqrt(2*np.pi*sigma1**2) * np.exp(-(x-ev1)**2/(2*sigma1**2)) + amp2/np.sqrt(2*np.pi*sigma2**2) * np.exp(-(x-ev2)**2/(2*sigma2**2))

def fg(datapointsX,datapointsY,amp, ev, sigma):
    assumtionsGauss = np.array([amp, ev, sigma])

    poptg, pcovg = curve_fit(gaussfitfunc, datapointsX, datapointsY,assumtionsGauss)

    return poptg
    
def fg2(datapointsX,datapointsY,amp1, ev1, sigma1, amp2, ev2, sigma2):
    assumtionsGauss = np.array([amp1, ev1, sigma1, amp2, ev2, sigma2])
    
    poptg, pcovg = curve_fit(gaussfitfunc2, datapointsX, datapointsY,assumtionsGauss)

    return poptg

def cl(chosenvars):
    varlist = list(globals().keys()) #all variable names are given in the dictionary globals which can be made a list by first use .keys() and then converted by list()
    returnvars = []
    if chosenvars[-1] == "*":
        varsidentifier = chosenvars[0:-1]
        for Numdata in range(len(varlist)):
            if varlist[Numdata][0:len(varsidentifier)] == varsidentifier:
                exec("returnvars.append(%s)"%varlist[Numdata])
    return returnvars

def cx(lst):
    return np.linspace(min(lst),max(lst),len(lst))

def d0(delta0fit):
    returnd0 = []
    if type(delta0fit) == list:
        for i in range(len(delta0fit)):
            returnd0i=[]
            for ii in range(len(delta0fit[i])):
                returnd0i.append(delta0fit[i][ii]-.001*ii)
            returnd0.append(returnd0i)
    else:
        for ii in range(len(delta0fit)):
            returnd0.append(delta0fit[ii]-.001*i)
    return returnd0
        
def pl(ydata):
    figA = plt.figure("Amplitude delta Correlation")
    axA = figA.add_subplot(111)

    if type(ydata[0]) == list:
        for Numdata in range(len(ydata)):
            axA.plot(range(len(ydata[Numdata])),ydata[Numdata])
    else:
        axA.plot(range(len(ydata)),ydata, 'g')
    if dp:
        plt.show()
    else:
        print("dp is set to 0 use plt.show()")

def pf(xdata,func):
    figA = plt.figure("Amplitude delta Correlation")
    axA = figA.add_subplot(111)
    axA.plot(xdata,func)
    if dp:
        plt.show()
    else:
        print("dp is set to 0 use plt.show()")

def cp(ydata1,ydata2):
    figA = plt.figure("Amplitude delta Correlation")
    axA = figA.add_subplot(111)
    #for Numdata in range(len(ydata1)):
    colors = ["r","g","b"]

    if type(ydata1[0]) == list:
        for Numdata in range(len(ydata1)):
            axA.plot(ydata1[Numdata],ydata2[Numdata], colors[Numdata]+'.')
    else:
        axA.plot(ydata1,ydata2, 'g.')
    if dp:
        plt.show()
    else:
        print("dp is set to 0 use plt.show()")

def c2(x1,y1,x2,y2): # conjunction2D
    #Using your x and y
    c1 = np.array([x1,y1]).T
    print(c1)
    c2 = np.array([x2,y2]).T

    radius = 0.0000035
    r = radius
    min_d = (2*r)*(2*r)

    shiftx = x1-np.mean(x1)
    shifty = y1-np.mean(y1)

    factor = max(shiftx)/max(shifty)
    print(factor)

    c1 = np.array([x1/factor,y1]).T
    c2 = np.array([x2/factor,y2]).T

    d = distance.cdist(c1,c2,'sqeuclidean')

    c1 = np.array([x1,y1]).T
    c2 = np.array([x2,y2]).T

    intersect = d <= min_d

    a,b = np.where(intersect)

    print(np.mean(np.transpose(c1[a])[0]))
    print(np.mean(np.transpose(c1[a])[1]))
    print(np.mean(np.transpose(c2[b])[0]))
    print(np.mean(np.transpose(c2[b])[1]))
    
    return np.transpose(c1[a])[0],np.transpose(c1[a])[1],a,b

def ht(ydata,num_bins,fc):
    n, bins, patches = plt.hist(ydata, num_bins, facecolor = fc)# = 'green')
    if dp:
        plt.show()
    else:
        print("dp is set to 0 use plt.show()")
    return bins[0:-1], n

def ld(chosenfolder,DataPath):
    filelist = os.listdir(DataPath)
    folderlist = []
    #folderoffset = "Withmathposition//"
    folderoffset = ""
    for i in range(len(filelist)):
        if "." not in filelist[i]:
            folderlist.append(filelist[i])
    if chosenfolder[0] == "*":
        if chosenfolder[-1] == "*":
            listidentifier = chosenfolder[1:-1]


            for Numdata in range(len(chosenfolders)):
                for i in range(len(targetfiles)):
                    if targetfiles[i][0:-4] == listidentifier and len(targetfiles[i][0:-4]) == len(listidentifier):
                        chosenset = listidentifier+chosenfolders[Numdata]
                        exec("global %s ; %s=[]"%(chosenset,chosenset))

                        with open(DataPath + "/All_Modular/" + folderoffset + chosenfolders[Numdata] + "/" + targetfiles[i],"rt") as csvfile:
                            spamreader = csv.reader(csvfile, delimiter=',')
                            for row in spamreader:
                                content = np.array(row).astype(float)
                                exec("%s.append(%10.10f)"%(chosenset,content))
                    else:
                        continue
    else:
        if type(chosenfolder) != str:
            for Numdata in range(len(chosenfolder)):
                for i in range(len(targetfiles)):
                    lstname = targetfiles[i][0:-4]
                    chosenset = lstname+chosenfolder[Numdata]
                    exec("global %s ; %s=[]"%(chosenset,chosenset))

                    with open(DataPath + "//All_Modular//" + folderoffset + chosenfolders[Numdata] + "//" + targetfiles[i],"rt") as csvfile:
                        spamreader = csv.reader(csvfile, delimiter=',')
                        for row in spamreader:
                            content = np.array(row).astype(float)
                            exec("%s.append(%10.10f)"%(chosenset,content))
        else:
            for i in range(len(targetfiles)):
                lstname = targetfiles[i][0:-4]
                chosenset = lstname+chosenfolder
                exec("global %s ; %s=[]"%(chosenset,chosenset))
                with open(DataPath + "//" + chosenfolder + "//" + targetfiles[i],"rt") as csvfile:
                    spamreader = csv.reader(csvfile, delimiter=',')
                    for row in spamreader:
                        content = np.array(row).astype(float)
                        exec("%s.append(%10.10f)"%(lstname+chosenfolder,content))
                        

def rmse(meangamma,gammalst):
    return np.sqrt(sum((meangamma-gammalst)**2.)/len(gammalst))

def means(choice):
    meanvalues = []
    meanvalueserr = []
    for Numdata in range(len(chosenfolders)):
        for i in range(len(targetfiles)):
            lstname = targetfiles[i][0:-4]
            chosenset = lstname+chosenfolders[Numdata] + str(Numdata)
            exec("%s=[]"%(chosenset))
            with open(DataPath + "/All_Modular/" + folderoffset + chosenfolders[Numdata] + "//" + targetfiles[i],"rt") as csvfile:
                spamreader = csv.reader(csvfile, delimiter=',')
                for row in spamreader:
                    content = np.array(row).astype(float)
                    exec("%s.append(%10.10f)"%(chosenset,content))
                if lstname[-6:0] != "stderr":
                    print(lstname)

    figA = plt.figure("meanvalues and errorbars")
    axA = figA.add_subplot(111)
    plt.errorbar(range(len(meanvalues)), meanvalues, xerr = np.array(meanvalueserr)*0., yerr = meanvalueserr) 
    plt.show()

def between(listtotest,lower,higher):
    listbetween = []
    for i in range(len(listtotest)):
        if listtotest[i] > lower and listtotest[i] < higher:
            listbetween.append(listtotest[i])
    return listbetween