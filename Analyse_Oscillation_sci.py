################################################################################################################
# Name: AnalysePhaseOnePixelLine.py
#
# Purpose: At the moment the most generalized version of the evaluationscript
#			
# Things to fix: Comments
#
#################################################################################################################

import matplotlib.pyplot as plt
import numpy as np
from numpy import transpose as nptp
import csv
import os, sys
import UndulatorDiffractionOneDimensional as UnDi
from natsort import natsorted, ns

from scipy.optimize import curve_fit, minimize_scalar

arg = sys.argv

def sinefit(x,a,lambdaOsz,phi,b,c):
	return (a * np.sin(2*np.pi*x/lambdaOsz + phi) + b ) * (1-c * x)

def sinefitgauss(x,a,lambdaOsz,phi,b,c,d):
	return (a * np.sin(2*np.pi*x/lambdaOsz + phi) + b ) * np.exp(-(x)**2/11.)
    
DataPath = arg[1]
print(DataPath)

#chosenfolder = "SharedPinMM"
chosenfolder = "SharedP"

AmplfitLst=[]
with open(DataPath + "//" + chosenfolder + "//" + "Amplfit" + ".csv","rt") as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    for row in spamreader:
        AmplfitLst.append(np.array(row).astype(np.float))

YtranLst=[]
with open(DataPath + "//" + chosenfolder + "//" + "Ytran" + ".csv","rt") as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    for row in spamreader:
        YtranLst.append(np.array(row).astype(np.float))

deltadLst=[]
with open(DataPath + "//" + chosenfolder + "//" + "deltad" + ".csv","rt") as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    for row in spamreader:
        deltadLst.append(np.array(row).astype(np.float))


AmplfitYtranLst=[]
with open(DataPath + "//" + chosenfolder + "//" + "AmplfitYtran" + ".csv","rt") as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    for row in spamreader:
        AmplfitYtranLst.append(np.array(row).astype(np.float))
        
uncertAmplLst=[]
with open(DataPath + "//" + chosenfolder + "//" + "uncertAmpl" + ".csv","rt") as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    for row in spamreader:
        uncertAmplLst.append(np.array(row).astype(np.float))
        
uncertYtranLst=[]
with open(DataPath + "//" + chosenfolder + "//" + "uncertYtra" + ".csv","rt") as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    for row in spamreader:
        uncertYtranLst.append(np.array(row).astype(np.float))
"""
AmplfitdeltadLst=[]
with open(DataPath + "//" + chosenfolder + "//" + "Amplfitdeltad" + ".csv","rt") as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    for row in spamreader:
        AmplfitdeltadLst.append(np.array(row).astype(np.float))

YtrandeltadLst=[]
with open(DataPath + "//" + chosenfolder + "//" + "Ytrandeltad" + ".csv","rt") as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    for row in spamreader:
        YtrandeltadLst.append(np.array(row).astype(np.float))
"""
"""
figA = plt.figure("Amplitude delta Correlation")
axA = figA.add_subplot(111)
axA.plot(range(len(AmplfitLst)),AmplfitLst, 'r')

figB = plt.figure("Amplitude offset Correlation")
axB = figB.add_subplot(111)
axB.plot(range(len(YtranLst)),YtranLst, 'g')

figC = plt.figure("sFactor delta Correlation")
axC = figC.add_subplot(111)
axC.plot(range(len(deltadLst)),deltadLst, 'b')

#figC = plt.figure("Comparison of Correlations")
#axC = figC.add_subplot(111)
#axC.plot(range(len(AmpldeltaCorrelLst)),AmpldeltaCorrelLst, 'r')
#axC.plot(range(len(AmploffsetCorrelLst)),AmploffsetCorrelLst, 'g')
#axC.plot(range(len(sFactordeltaLst)),sFactordeltaLst, 'b')
plt.show()
#"""
yy = 0.
#print(nptp(YtranLst)[0])


AmplfitLst=nptp(AmplfitLst)
YtranLst=nptp(YtranLst)
deltadLst=nptp(deltadLst)

AmplfitYtranLst=nptp(AmplfitYtranLst)
#AmplfitdeltadLst=nptp(AmplfitdeltadLst)
#YtrandeltadLst=nptp(YtrandeltadLst)

print(AmplfitYtranLst)

uncertAmplLst=nptp(uncertAmplLst)
uncertYtranLst=nptp(uncertYtranLst)

figV = plt.figure("Amplitun")
axA = figV.add_subplot(111)
axA.plot(range(len(uncertAmplLst[0])),np.array(uncertAmplLst[0]), 'r')
axA.plot(range(len(uncertYtranLst[0])),np.array(uncertYtranLst[0])*1000., 'g')
axA.plot(range(len(uncertYtranLst[0])),np.array(uncertAmplLst[0])*np.array(uncertYtranLst[0])*10000., 'g')
#axA.plot(range(len(lowcorrelYtranLst)),lowcorrelYtranLst, 'g')
#axA.plot(range(len(lowcorreldeltadLst)),lowcorreldeltadLst, 'b')

plt.show()

figV = plt.figure("Amplitun")
axA = figV.add_subplot(111)
axA.plot(range(len(YtranLst[0])),np.array(YtranLst[0])*1000., 'g')
#axA.plot(range(len(lowcorrelYtranLst)),lowcorrelYtranLst, 'g')
#axA.plot(range(len(lowcorreldeltadLst)),lowcorreldeltadLst, 'b')

#plt.show()

figV = plt.figure("Amplitun")
axA = figV.add_subplot(111)
axA.plot(range(len(AmplfitYtranLst[0])),AmplfitYtranLst[0], 'r')
#axA.plot(range(len(AmplfitdeltadLst[0])),AmplfitdeltadLst[0], 'g')
#axA.plot(range(len(YtrandeltadLst[0])),YtrandeltadLst[0], 'b')

#axA.plot(range(len(AmplfitdeltadLst[0])),AmplfitYtranLst[0]*AmplfitdeltadLst[0]*YtrandeltadLst[0], 'k')

lowcorrelAmplfitLst = []
lowcorrelYtranLst = []
lowcorreldeltadLst = []

for i in range(len(AmplfitLst[0])):
    if abs(np.array(uncertAmplLst[0][i])*np.array(uncertYtranLst[0][i])) < 10./10000000.:#if AmplfitYtranLst[0][i] > -.001 and AmplfitYtranLst[0][i] < .0:
        if abs(AmplfitYtranLst[0][i]) < .01:#if AmplfitYtranLst[0][i] > -.001 and AmplfitYtranLst[0][i] < .0:
            lowcorrelAmplfitLst.append(AmplfitLst[0][i])
            lowcorrelYtranLst.append(YtranLst[0][i])
            print("This value: ",i)
            #lowcorreldeltadLst.append(deltadLst[0][i])


plt.show()

figV = plt.figure("Amplitun")
axA = figV.add_subplot(111)
axA.plot(range(len(lowcorrelAmplfitLst)),lowcorrelAmplfitLst, 'r')
axA.plot(range(len(lowcorrelYtranLst)),np.array(lowcorrelYtranLst)*1000., 'g')
#axA.plot(range(len(lowcorreldeltadLst)),lowcorreldeltadLst, 'b')

plt.show()


print(YtranLst)

oscIntensity = []

z0 = 10.
#"""
for i in range(len(AmplfitLst[0])):
    oscIntensity.append(UnDi.diffractionintensity([yy]*len(YtranLst[0]),YtranLst[0][i],0,AmplfitLst[0][i],z0,404.3e-9,deltadLst[0][i]).real)

figA = plt.figure("Amplitude delta Correlation")
axA = figA.add_subplot(111)
axA.plot(range(len(oscIntensity)),oscIntensity, 'r')
#print(oscIntensity)
plt.show()
#"""

deltadLst0=[]

for i in range(len(AmplfitLst[0])):
    deltadLst0.append((deltadLst[0][i]-.245)/((i+1)*.001))
    print(deltadLst[0][i])
    

n, bins, patches = plt.hist(deltadLst0, 2500, facecolor='g')


plt.xlabel('Smarts')
plt.ylabel('Probability')
#plt.title('Histogram of IQ')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
#plt.xlim(40, 160)
#plt.ylim(0, 0.03)
plt.grid(True)
plt.show()

"""
rat4 = 1.65
a2 = .004
NumOfstripes = 5
lambdaU = .08
zetaprime2=np.linspace( -(rat4)*a2,(rat4)*a2,1000)

print("Now Undupos")
print(len(UnduposLst))

print("Now amplitude")
print(len(amplitudeLst))

#print("Now sFactorLst")
#print(len(sFactorLst))

print("Now asymmetryLst")
print(len(asymmetryLst))

print("Now offsets")
print(len(offsetLst))

print("Now wavelengths")
print(len(center0Lst))

OscLst = []

#center0 = 4.088331686445741e-07

SourcelstIntensityAllPos1=[]
SourcelstIntensityAllPos2=[]
SourcelstIntensityAllPos3=[]
SourcelstIntensityAllPos4=[]
SourcelstIntensityAllPos5=[]

sFactor = 9.5

for Stripeindex in range(NumOfstripes):
    OscLstStripe=[]
    OscLstStripe1=[]

    for Undulatorpositionindex in range(0,int(825/1)):
        Uix = Undulatorpositionindex

        WLofPixel = 409.140e-9 +Stripeindex*-.478e-9
        center0 = center0Lst[Uix]
        print(int((UnduposLst[Uix][Stripeindex]-0.)/5))
        print(int(Uix))
        if 1:#int((UnduposLst[Uix][Stripeindex]-0.)/5) == int(Uix):####20.

            Sourcelst=[]
            for iii in range(0,int(len(zetaprime2)/1),10):
                lambdaNull = Undu.lambda0fit(center0,(zetaprime2[iii]-offsetLst[int(Uix)][0])/sFactor,lambdaU)
                deltaL = WLofPixel-lambdaNull
                if (zetaprime2[iii]-asymmetryLst[Uix] > -a2 and zetaprime2[iii]-asymmetryLst[Uix] < a2):
                    Sourcelst.append(amplitudeLst[int(Uix)][0]*Undu.analyticalspecPhase((zetaprime2[iii]-offsetLst[int(Uix)][0])/sFactor,deltaL,center0,WLofPixel,deltadLst[int(Uix)][0]))
                    #if zetaprime2[iii]==zetaprime2[500]:
                        #print(Sourcelst[-1])
                        #print(Uix)
                else:
                    Sourcelst.append(0.)

        else:
            if int(UnduposLst[Uix][Stripeindex]) == int(-1):
                Sourcelst=[]
                for iii in range(0,len(zetaprime2),1):
                    Sourcelst.append(-1)

            Sourcelst = nptp(Sourcelst)

        SourcelstIntensity = abs(np.array(Sourcelst)**2.)
        if Stripeindex==0:
            SourcelstIntensityAllPos1.append(SourcelstIntensity)
        if Stripeindex==1:
            SourcelstIntensityAllPos2.append(SourcelstIntensity)
        if Stripeindex==2:
            SourcelstIntensityAllPos3.append(SourcelstIntensity)
        if Stripeindex==3:
            SourcelstIntensityAllPos4.append(SourcelstIntensity)
        if Stripeindex==4:
            SourcelstIntensityAllPos5.append(SourcelstIntensity)

        if Undulatorpositionindex+2 > len(UnduposLst):
            print("exit!")
            break

print("now plot")
posplt = False

if posplt:
    for i in range(len(SourcelstIntensityAllPos1)):
        figA = plt.figure("Oscillation")
        axA = figA.add_subplot(111)
        if UnduposLst[i][0] != -1:
            axA.plot(zetaprime2,SourcelstIntensityAllPos1[i])
        if UnduposLst[i][1] != -1:
            axA.plot(zetaprime2,SourcelstIntensityAllPos2[i])
        if UnduposLst[i][2] != -1:
            axA.plot(zetaprime2,SourcelstIntensityAllPos3[i])
        if UnduposLst[i][3] != -1:
            axA.plot(zetaprime2,SourcelstIntensityAllPos4[i])
        if UnduposLst[i][4] != -1:
            axA.plot(zetaprime2,SourcelstIntensityAllPos5[i])
        plt.show()

print("now plot")
oscplot1=[]
oscplotX1=[]
oscplot2=[]
oscplotX2=[]
oscplot3=[]
oscplotX3=[]
oscplot4=[]
oscplotX4=[]
oscplot5=[]
oscplotX5=[]
figB = plt.figure("Oscillation")
axB = figB.add_subplot(111)
for i in range(len(SourcelstIntensityAllPos1)):
    if UnduposLst[i][0] != -1:
        oscplotX1.append(i)
        oscplot1.append(SourcelstIntensityAllPos1[i][50])
    if UnduposLst[i][1] != -1:
        oscplotX2.append(i)
        oscplot2.append(SourcelstIntensityAllPos2[i][50])
    if UnduposLst[i][2] != -1 and SourcelstIntensityAllPos3[i][50] < .004: # .003:
        oscplotX3.append(UnduposLst[i][0]/1000.)

        oscplot3.append(SourcelstIntensityAllPos3[i][50])
    if UnduposLst[i][3] != -1:
        oscplotX4.append(i)
        oscplot4.append(SourcelstIntensityAllPos4[i][50])
    if UnduposLst[i][4] != -1:
        oscplotX5.append(i)
        oscplot5.append(SourcelstIntensityAllPos5[i][50])


"""
"""
axB.plot(oscplotX1,oscplot1, 'r.')
axB.plot(oscplotX2,oscplot2, 'g.')
axB.plot(oscplotX3,oscplot3, 'b.')
axB.plot(oscplotX4,oscplot4, 'c.')
axB.plot(oscplotX5,oscplot5, 'k.')
plt.show()
"""
"""
#axB.plot(oscplotX1,oscplot1, 'r')
#axB.plot(oscplotX2,oscplot2, 'g')
axB.plot(oscplotX3,oscplot3, 'b')
#axB.plot(oscplotX4,oscplot4, 'c')
#axB.plot(oscplotX5,oscplot5, 'k')
plt.show()
"""
