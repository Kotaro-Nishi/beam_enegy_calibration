################################################################################################################
# Name: Analyse_gamma.py
#
# Purpose: Plot Fit results
#			
# Things to fix: Comments
#
#################################################################################################################

import matplotlib.pyplot as plt
import numpy as np
from numpy import transpose as nptp
import csv
import os, sys
from natsort import natsorted, ns
from scipy.optimize import curve_fit
import pandas as pd

from scipy.spatial import distance

arg = sys.argv

import random
import statistics

global dp
dp = 1 # do plot instant

"""
# Creating bins
x_min = np.min(x)
x_max = np.max(x)
  
y_min = np.min(y)
y_max = np.max(y)
  
x_bins = np.linspace(x_min, x_max, 50)
y_bins = np.linspace(y_min, y_max, 20)
  
fig, ax = plt.subplots(figsize =(10, 7))
# Creating plot
plt.hist2d(x, y, bins =[x_bins, y_bins])
plt.title("Changing the bin scale")
  
ax.set_xlabel('X-axis') 
ax.set_ylabel('X-axis') 
  
# show plot
plt.tight_layout() 
plot.show()
#"""

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
            #print(varlist[Numdata][0:len(varsidentifier)])
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

    if type(ydata1[0]) == list:
        for Numdata in range(len(ydata1)):
            axA.plot(ydata1[Numdata],ydata2[Numdata], '.')
    else:
        axA.plot(ydata1,ydata2, 'b.')
    plt.show()

def c2(x1,y1,x2,y2): # conjunction2D
    #Using your x and y
    c1 = np.array([x1,y1]).T
    c2 = np.array([x2,y2]).T

    radius = 0.0000054
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
    
    return np.transpose(c1[a])[0],np.transpose(c1[a])[1],a,b#,np.transpose(c2[b])[0],np.transpose(c2[b])[1]

def ht(ydata,num_bins):
    n, bins, patches = plt.hist(ydata, num_bins, facecolor = 'green')
    if dp:
        plt.show()
    else:
        print("dp is set to 0 use plt.show()")
    return bins[0:-1], n
        
def ld(chosenfolder):
    filelist = os.listdir(DataPath)
    folderlist = []
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

                        with open(DataPath + "//" + chosenfolders[Numdata] + "//" + targetfiles[i],"rt") as csvfile:
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

                    with open(DataPath + "//" + chosenfolder[Numdata] + "//" + targetfiles[i],"rt") as csvfile:
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

chosenfolders = [
                #"SharedP_20210401_100_000425_825_seed10000_new_zpfit",
                ###"SharedP_20210401_100_000425_825_seed10000_new_zpfit_NoConv460",
                "SharedP_20210401_100_000425_825_seed10000_new_zpfit_NoConv",
                #"SharedP_20210401_100_025450_825_seed10000_new_zpfit",
                #"SharedP_20210401_100_050475_825_seed10000_new_zpfit",
                #"SharedP_20210401_100_075500_825_seed10000_new_zpfit",
                #"SharedP_20210401_100_100525_825_seed10000_new_zpfit",
                #"SharedP_20210401_100_125550_825_seed10000_new_zpfit",
                #"SharedP_20210401_100_150575_825_seed10000_new_zpfit",
                #"SharedP_20210401_100_175600_825_seed10000_new_zpfit",
                #"SharedP_20210401_100_200625_825_seed10000_new_zpfit",
                #"SharedP_20210401_100_225650_825_seed10000_new_zpfit",
                #"SharedP_20210401_100_250675_825_seed10000_new_zpfit",
                #"SharedP_20210401_100_275700_825_seed10000_new_zpfit",
                #"SharedP_20210401_100_300725_825_seed10000_new_zpfit",
                #"SharedP_20210401_100_325750_825_seed10000_new_zpfit",
                #"SharedP_20210401_100_350775_825_seed10000_new_zpfit",
                #"SharedP_20210401_100_375800_825_seed10000_new_zpfit",
                
                #"SharedP_20210401_100_400825_825_seed10000_new_zpfit",
                
                ###"SharedP_20210401_100_400825_825_seed10000_new_zpfit_NoConv460",
                
                "SharedP_20210401_100_400825_825_seed10000_new_zpfit_NoConv",                
                
                #"SharedP_20210401_100_000275_825_seed10000_new_zpfit",
                #"SharedP_20210401_100_275550_825_seed10000_new_zpfit",
                #"SharedP_20210401_100_550825_825_seed10000_new_zpfit",
                #"SharedP_20210401_100_000275_825_seed10000_new_zpfit",
                #"SharedP_20210401_825_825_seed10000_new_zpfit_x",
                ###"SharedP_20210401_825_825_seed10000_new_zpfit_x_NoConv460",
                "SharedP_20210401_825_825_seed10000_new_zpfit_x_NoConv",
                #"SharedP_20210401_825_seed10000_MultiWL",
                "SharedP_20210401_825_seed10000_MultiWL_NoConv",

                "SharedP_20210401_825_825_seed10000_new_zpfit_x_NoConv",
                ]
                
                


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
            #exec("print(%s%d)"%(chosenset,Numdata))
            with open(DataPath + "//" + chosenfolders[Numdata] + "//" + targetfiles[i],"rt") as csvfile:
                spamreader = csv.reader(csvfile, delimiter=',')
                for row in spamreader:
                    content = np.array(row).astype(float)
                    exec("%s.append(%10.10f)"%(chosenset,content))
            #print("gfdgfsdg" + lstname[-6])
                if lstname[-6:0] != "stderr":#"gammafitSharedP":
                    print(lstname)
                    #exec("print(np.mean(%s))"%np.array(chosenset))
                    #exec("print(rmse(np.mean(%s),%s))"%(np.array(chosenset),np.array(chosenset)))
                    #exec("meanvalues.append(np.mean(%s))"%np.array(chosenset))
                    #exec("meanvalueserr.append(rmse(np.mean(%s),%s))"%(np.array(chosenset),np.array(chosenset)))

    #pl(meanvalues)
    figA = plt.figure("meanvalues and errorbars")
    axA = figA.add_subplot(111)
    plt.errorbar(range(len(meanvalues)), meanvalues, xerr = np.array(meanvalueserr)*0., yerr = meanvalueserr) 
    plt.show()
                #exec("print(%s)"%chosenset)

def between(listtotest,lower,higher):
    listbetween = []
    for i in range(len(listtotest)):
        if listtotest[i] > lower and listtotest[i] < higher:
            listbetween.append(listtotest[i])
    return listbetween
    




random.seed(10000)
NumOfSeries = 825




"""
40
66
69
77
88
119
147
196
"""

"""

for ii in range(300):
    listOfboots = []
    for i in range(NumOfSeries):
        listOfboots.append(random.randint(0,NumOfSeries-1))
    print(np.average(np.array(listOfboots),weights=np.array(listOfboots)))
    print(np.average(np.array(listOfboots)))
    print(min(listOfboots))
    print(max(listOfboots))
    print(statistics.median(listOfboots))
    if ii > 117:
        ht(listOfboots,825)
#"""
DataPath = arg[1]
print(DataPath)

chosenfolder = "Interval700"#"Interval700"

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
               "zpfit.csv"]
"""
for i in range(len(targetfiles)):
    lstname = targetfiles[i][0:-4]
    exec("%s=[]"%(lstname+chosenfolder))

    with open(DataPath + "//" + chosenfolder + "//" + targetfiles[i],"rt") as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        for row in spamreader:
            content = np.array(row).astype(float)
            exec("%s.append(%10.10f)"%(lstname+chosenfolder,content))
            
"""
figA = plt.figure("Amplitude delta Correlation")
axA = figA.add_subplot(111)
#axA.plot(range(len(Amplfit)),Amplfit, 'r')
#axA.plot(range(len(delta0fit)),Ytran, 'r')
#axA.plot(range(len(bfit)),delta0fit, 'r')
#axA.plot(range(len(bfit)),bfit, 'r')
axA.plot(range(len(bfit)),gammafit, 'g')
#axA.plot(range(len(bfit)),bfit, 'r')
#axA.plot(range(len(bfit)),bfit, 'r')

#"""
chosenfolder = "Interval725"#"SharedP"#_10_with_delta"#_700_11m"

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
               "delta0fitcfitCorrel.csv"]

for i in range(len(targetfiles)):
    lstname = targetfiles[i][0:-4]
    exec("%s=[]"%(lstname))

    with open(DataPath + "//" + chosenfolder + "//" + targetfiles[i],"rt") as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        for row in spamreader:
            content = np.array(row).astype(float)
            exec("%s.append(%10.10f)"%(lstname,content))
#"""
figA = plt.figure("Amplitude delta Correlation")
#axA = figA.add_subplot(111)
axA.plot(range(len(gammafit)),np.array(gammafit), 'b')
#axA.plot(range(len(Ytran)),Ytran, 'r')
#axA.plot(range(len(gammafitstderr)),gammafitstderr, 'r')

plt.show()

#figB = plt.figure("Amplitude offset Correlation")
#axB = figB.add_subplot(111)
#axB.plot(range(len(AmploffsetCorrelLst)),AmploffsetCorrelLst, 'g')

cg04,cd04,aa,bb = d2(gn0[0],dn0[0],gn4[0],dn4[0])
cp([gn0[0],gn4[0],gn8[0],[352.414],cg04],[dn0[0],dn4[0],dn8[0],[.1870],cd04])

ld("*delta0fit*")
ld("*gammafit*")

gn0 = cl("gammafitSharedP_20210401_100_000425_825_seed10000_new_zpfit_NoConv*")
gn4 = cl("gammafitSharedP_20210401_100_400825_825_seed10000_new_zpfit_NoConv*")
dn0 = cl("delta0fitSharedP_20210401_100_000425_825_seed10000_new_zpfit_NoConv*")
dn4 = cl("delta0fitSharedP_20210401_100_400825_825_seed10000_new_zpfit_NoConv*")
gn8 = cl("gammafitSharedP_20210401_825_825_seed10000_new_zpfit_x_NoConv*")
dn8 = cl("delta0fitSharedP_20210401_825_825_seed10000_new_zpfit_x_NoConv*")



gg0 = cl("gammafitSharedP_20210401_100_000425_825_seed10000_new_zpfit*")
gg4 = cl("gammafitSharedP_20210401_100_400825_825_seed10000_new_zpfit*")
gg8 = cl("gammafitSharedP_20210401_825_825_seed10000_new_zpfit_x*")
dd0 = cl("delta0fitSharedP_20210401_100_000425_825_seed10000_new_zpfit*")
dd4 = cl("delta0fitSharedP_20210401_100_400825_825_seed10000_new_zpfit*")
dd8 = cl("delta0fitSharedP_20210401_825_825_seed10000_new_zpfit_x*")
gwl = cl("gammafitSharedP_20210401_825_seed10000_MultiWL_NoConv*")

gn0 = cl("gammafitSharedP_20210401_100_000425_825_seed10000_new_zpfit_NoConv460*")
gn4 = cl("gammafitSharedP_20210401_100_400825_825_seed10000_new_zpfit_NoConv460*")
dn0 = cl("delta0fitSharedP_20210401_100_000425_825_seed10000_new_zpfit_NoConv460*")
dn4 = cl("delta0fitSharedP_20210401_100_400825_825_seed10000_new_zpfit_NoConv460*")
gn8 = cl("gammafitSharedP_20210401_825_825_seed10000_new_zpfit_x_NoConv460*")
dn8 = cl("delta0fitSharedP_20210401_825_825_seed10000_new_zpfit_x_NoConv460*")