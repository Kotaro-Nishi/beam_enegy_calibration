import ModularFunctions as MoFu

import matplotlib.pyplot as plt
import numpy as np
from numpy import transpose as nptp
import csv
import os, sys
from natsort import natsorted, ns
from scipy.optimize import curve_fit

from scipy.spatial import distance

folder= "D:\Doktor\BeamtimeFeb2020_Common\data\Beam20200320_1\HeidenheinBeam20200320_1"

filename = "pos0.txt"

print(folder)

filelist = os.listdir(folder)

filelist = natsorted(filelist)

print(filelist)

poslist=[]

#works for positions smaller 1000.00000 mm
for f in range(len(filelist)):
    with open(folder+ "//" + filelist[f], "r") as csv_file:
        csvdata = list(csv.reader(csv_file, delimiter=','))
        csvdata = csvdata[0]

        backnumber = []
        frontnumber = []
        precomma = 1
        
        for i in range(len(csvdata)):
            if csvdata[i]==" ":
                continue

            if csvdata[i]=="+":
                sign=1

            if csvdata[i]=="-":
                sign=-1

            if csvdata[i]==".":
                precomma=0

            if csvdata[i].isdigit() and precomma==1:
                frontnumber.append(csvdata[i])

            if csvdata[i].isdigit() and precomma==0:
                backnumber.append(csvdata[i])

        number = "".join(frontnumber+backnumber)
        #print(sign*float(number)/100000.)
        pos = sign*float(number)/100000.#+f

    poslist.append(pos)

fntsze = 2*13
figA = plt.figure("Positions")
axA = figA.add_subplot(111)


axA.plot(list(range(len(poslist))), poslist,'r')

axA.set_xlabel(r'Undulator position $\Delta$ d [cm]',fontsize=fntsze)
axA.set_ylabel('Intensity [a.u.]',fontsize=fntsze)
axA.tick_params(axis='both', which='major', labelsize=fntsze)
axA.legend(['Fit','Data'],fontsize=fntsze)
plt.show()

print(poslist)

"""
with open(DataPath + SaveSharedP + "\\zpfitelevationfitCorrel.csv", 'a+', newline='') as write_obj:
	# Create a writer object from csv module
	csv_writer = csv.writer(write_obj)
	# Add contents of list as last row in the csv file
	csv_writer.writerow([zpfitelevationfitCorrel])
"""