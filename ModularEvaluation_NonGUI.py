import matplotlib.pyplot as plt
import UndulatorDiffractionOneDimensional_correct_XAmpl_InitialFit_NoConvolve as UnDi
import InterferenceFunction_more_x as IFunc
import AmplitudeFormula_x as AUndu
import importlib

import os

#open("./saving_things.py").read()

"""
try:
    import saving_things
except:
"""
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

ratio = 1000. / 606#606. #how many pixel are aperture
UnDi.apy = .0081/2.
zetaxx = np.array([0.])

zetaprime2 = np.linspace(-ratio*UnDi.apy,ratio*UnDi.apy,1000)

NumOfSeries = 825
LengthOfData = 825
Interval = [000,825]
#Interval = [000,425]

WidthOfStripe = 5
print(WidthOfStripe)

dispersion = 0.004808759252942803e-9 #Former 0.0048085643228767475e-9
pixelof404nm = 1315.3024745554976 #Former: 1243.8551272671261
pixelof404nm = pixelof404nm/WidthOfStripe
CalibWl = 404.6565e-9

#(407.7837-404.6565)

#global listOfboots
#global pltweights
#listOfboots = []
#pltweights = []

Yscreenfit = (330)*10**-6
LambdaRfit = 412e-9

def per_iteration(pars, iteration, resid, *args, **kws):
    
    print(" resid ", sum(resid))
    print(" residSqr ", sum(resid**2.))

    print("Amplfit: " + str(pars['Amplfit_%i'% listOfboots[0]].value))
    print("zpfit: " + str(pars['zpfit_%i'% listOfboots[0]].value))
    print('Ytran: ' + str(pars['Ytran_%i'% listOfboots[0]].value))
    print("gammafit: " + str(pars['gamma_%i'% listOfboots[0]].value))
    print("delta0fit: " + str(pars['delta0_%i'% listOfboots[0]].value))
    print("elevationfit: " + str(pars['elevationfit_%i'% listOfboots[0]].value))
    print("cfit: " + str(pars['cfit_%i'% listOfboots[0]].value))

    Amplfit = pars['Amplfit_%i'% listOfboots[0]].value
    Ytran = pars['Ytran_%i'% listOfboots[0]].value
    gammafit = pars['gamma_%i'% listOfboots[0]].value
    delta0fit = pars['delta0_%i'% listOfboots[0]].value
    zpfit = pars['zpfit_%i'% listOfboots[0]].value
    elevationfit = pars['elevationfit_%i'% listOfboots[0]].value
    cfit = pars['cfit_%i'% listOfboots[0]].value

    UnDi.SingleAmplitudeDone = 0
    UnDi.AmplXY = []
    

def diffraction_dataset(params,i,yy):
    Amplfit = params['Amplfit_%i' % (i)].value
    Ytran   = params['Ytran_%i' % (i)].value
    gammafit  = params['gamma_%i' % (i)].value
    delta0fit  = params['delta0_%i' % (i)].value
    zpfit  = params['zpfit_%i' % (i)].value
    elevationfit  = params['elevationfit_%i' % (i)].value
    cfit = params['cfit_%i' % (i)].value
    bfit = 0
    Diff_intensity = UnDi.diffractionintensityInitialFit_X(yy,zetaxx,Ytran,i,AUndu.z0,WLofPixel,LambdaRfit,delta0fit+i*.001,gammafit,Amplfit,bfit,zpfit,Yscreenfit,elevationfit,whichdata,cfit)

    return Diff_intensity

def objective_diffraction(params, yy, data,weights):
    ndata,nx = data.shape
    
    resid = 0.0*data[:]

    for Strpidx in range(ndata):
        if Strpidx > 0:
            UnDi.SingleAmplitudeDone = 1

        diffractionresult = diffraction_dataset(params,listOfboots[Strpidx],yy)
        #1.5 is std of camera structure
        #1/.2 is conversion factor
        #1/2. is correction of convolutionsigma???Welches???
        #.05 makes strict positive
        #1*.2 is conversion factor
        #1/.7 is QE
        
        #1/5 average over five pixel columns
        #1/4 average over four images

        #resid[Strpidx] = (diffractionresult - data[Strpidx])*(1./(np.array([1.5/1.]*1000)+(1./.2)*(1/5.)*(1/4.)*np.sqrt(0.05+abs(diffractionresult)*.2/.7)))
        resid[Strpidx] = (diffractionresult - data[Strpidx])*(1./(np.array([3.5]*1000)+(np.sqrt(1/5.)*np.sqrt(1/4.)*(.7/.2)*np.sqrt(0.05+abs(diffractionresult)*.2/.7))))
    return np.array(resid.flatten())

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

whichdata = 0

print(whichdata)

slashes = [i for i,ltr in enumerate(DataPath) if ltr == "\\"]
DataPathBase = DataPath[0:slashes[-2]]

print(DataPathBase)

#datarangelow = 1804 + 550

Backgroundimage = cv2.imread(DataPathBase + "/BG/raw_shortExp/File1.tif",cv2.IMREAD_UNCHANGED)
ListOfAllMatrices = []
print(DataPathBase + "/BG/raw_shortExp/File1.tif")



List825VolumeMinusBG4Mean=[]

k0=0

for ii in range(k0,k0+LengthOfData):
    List4VolumeMinusBG = []
    print(filelist[ii*4])
    for i in range(4):

        PictureMinusBG = cv2.imread(DataPath + "/" + filelist[i+ii*4],cv2.IMREAD_UNCHANGED).astype(int)-Backgroundimage

        PictureMinusBG = nptp(nptp(np.float32(PictureMinusBG))[0:2300])
        
        PictureMinusBG = cv2.resize(PictureMinusBG,(int(2300/WidthOfStripe),1000),interpolation=cv2.INTER_AREA) #460*5=2300

        List4VolumeMinusBG.append(PictureMinusBG)

    PictureMinusBG4Mean = sum(List4VolumeMinusBG)/4.

    List825VolumeMinusBG4Mean.append(PictureMinusBG4Mean)

#SaveSharedP = "\All_Modular\SharedP_20210401_825_825_seed10000_new_zpfit_x"
#SaveSharedP = "\All_Modular\SharedP_20210401_825_825_seed10000_new_zpfit_boot000"
#SaveSharedP = "\All_Modular\SharedP_20210401_825_825_seed10000_new_zpfit_boot400"
#SaveSharedP = "\All_Modular\SharedP_20210401_825_825_seed10000_new_zpfit_boot800"
SaveSharedP = "\All_Modular\SharedP_20210401_825_825_seed10000_new_zpfit_MultiWL"

#D:\Doktor\BeamtimeFeb2020_Common\data\Beam20200320_1\Energymeasurement_2\raw_shortExp\All_Modular

print(NumOfSeries)
checkstart = 1
DoNotSave = 0

random.seed(10000)
#for i in range(NumOfSeries*6500):
#    random.randint(0,LengthOfData-1)

def ModFit(
zetaprime2,
List825VolumeMinusBG4Mean,
WidthOfStripe,
pixelof404nm,
dispersion,
CalibWl,
NumOfSeries,
LengthOfData,
Interval,
WLidxList,
PositionidxList):

    residsqrlst=[]
    print(len(WLidxList))
    if len(WLidxList) > 1 and len(PositionidxList) > 1:
        print("Wrong index argument")
        exit()
    else:
        if len(WLidxList) > 1:
            LoopidxList = WLidxList
            print("Wavelength variation")
        else:
            LoopidxList = PositionidxList
            print("Position variation")


    for Loopidx in LoopidxList:#range(LoopidxList[0],LoopidxList[-1]):
        print("Loopidx " + str(Loopidx))

        if len(WLidxList) > 1:
            chosenpx = Loopidx
        else:
            chosenpx = WLidxList[0]

        #chosenpx = 300
        #chosenpx = int(chosenpx/WidthOfStripe)

        FirstWL = (0 - pixelof404nm)*dispersion*WidthOfStripe + CalibWl
        global WLofPixel
        WLofPixel = 0
        WLofPixel = (chosenpx - pixelof404nm)*dispersion*WidthOfStripe + CalibWl

        #LastWL = (2304 - pixelof404nm)*dispersion*WidthOfStripe + CalibWl
        LastWL = (460 - pixelof404nm)*dispersion*WidthOfStripe + CalibWl

        print("First wavelength: ", FirstWL)
        print("chosenpx", chosenpx)
        print("Wavelength at which the analysis is done: ", WLofPixel)
        print("Last wavelength: ", LastWL)

        dataBooted = []
        global listOfboots
        listOfboots = []
        UnDi.SingleAmplitudeDone = 0
        for i in range(NumOfSeries):
            listOfboots=list(range(NumOfSeries))
            #listOfboots.append(random.randint(Interval[0],Interval[-1]-1))

        for SIt in range(len(listOfboots)):
            if listOfboots.count(listOfboots[SIt]) == 1:
                FirstSingularitem = SIt#listOfboots[SIt]
                break

        residLst = []

        interpolatedpicture = []

        for ii in range(LengthOfData):

            stripe = List825VolumeMinusBG4Mean[ii]

            stripe = nptp(stripe)[chosenpx]

            stripejjj=[]
            for jjj in range(len(stripe)):
                stripejjj.append(stripe[jjj])

            interpolatedpicture.append(stripe)

        data = np.array(interpolatedpicture)
        
        #print("stripe",stripe)

        for Boot in range(NumOfSeries):
            dataBooted.append(data[listOfboots[Boot]])
        dataBooted = np.array(dataBooted)

        dataLst = []
        for Strpidx in range(NumOfSeries):
            dataLst.append(dataBooted[Strpidx][500])
        
        """
        fntsze = 13
        figA = plt.figure("Diffraction")
        axA = figA.add_subplot(111)
        axA.plot(listOfboots, dataLst, 'b.')
        axA.set_xlabel('Position on Aperture [mm]',fontsize=fntsze)
        axA.set_ylabel('Intensity [a.u.]',fontsize=fntsze)
        axA.tick_params(axis='both', which='major', labelsize=fntsze)
        plt.show()
        #"""
        
        #Guessvalueslmfit = [5.4,.0001,352.4,.1870,5.3,15.,1.]
        #if Loopidx == LoopidxList[0]:
        Guessvalueslmfit = [1.4,.0001,352.1,.1870,5.,11.,1.]#old 5.2
        #else:
        #    Guessvalueslmfit = [1.4,.0001,352.1,deltaGuess,5.,11.,1.]#old 5.2

        fit_params = Parameters()

        for iy,y, in enumerate(dataBooted):
            fit_params.add('Amplfit_%i' % (listOfboots[iy]),value=Guessvalueslmfit[0])#, max=2.5, min=1.23)#, min = Guessvalueslmfit[0]-1.1, max=Guessvalueslmfit[0]+1.1)
            fit_params.add('Ytran_%i' % (listOfboots[iy]),value=Guessvalueslmfit[1])#, min = Guessvalueslmfit[1]-.0001, max=Guessvalueslmfit[1]+.0001)
            fit_params.add('gamma_%i' % (listOfboots[iy]),value=Guessvalueslmfit[2])
            fit_params.add('delta0_%i' % (listOfboots[iy]),value=Guessvalueslmfit[3])#, min = Guessvalueslmfit[4]-.1, max=Guessvalueslmfit[4]+.1)
            fit_params.add('zpfit_%i' % (listOfboots[iy]),value=Guessvalueslmfit[4])#, min = Guessvalueslmfit[4]-.1, max=Guessvalueslmfit[4]+.1)
            fit_params.add('elevationfit_%i' % (listOfboots[iy]),value=Guessvalueslmfit[5])#, min = Guessvalueslmfit[4]-.1, max=Guessvalueslmfit[4]+.1)
            fit_params.add('cfit_%i' % (listOfboots[iy]),value=Guessvalueslmfit[6])#, min = Guessvalueslmfit[4]-.1, max=Guessvalueslmfit[4]+.1)

        listOfbootsnoFirst = [x for x in listOfboots if x != listOfboots[FirstSingularitem]]
        for iy in tuple(range(len(listOfbootsnoFirst))):#tuple(range(NumOfSeries)[1:]):
            fit_params['Amplfit_%i' % listOfbootsnoFirst[iy]].expr='Amplfit_%i'% listOfboots[FirstSingularitem]
            
        for iy in tuple(range(len(listOfbootsnoFirst))):#tuple(range(NumOfSeries)[1:]):
            fit_params['Ytran_%i' % listOfbootsnoFirst[iy]].expr='Ytran_%i'% listOfboots[FirstSingularitem]
            
        for iy in tuple(range(len(listOfbootsnoFirst))):#tuple(range(NumOfSeries)[1:]):
            fit_params['gamma_%i' % listOfbootsnoFirst[iy]].expr='gamma_%i'% listOfboots[FirstSingularitem]
            
        for iy in tuple(range(len(listOfbootsnoFirst))):#tuple(range(NumOfSeries)[1:]):
            fit_params['delta0_%i' % listOfbootsnoFirst[iy]].expr='delta0_%i'% listOfboots[FirstSingularitem]

        for iy in tuple(range(len(listOfbootsnoFirst))):#tuple(range(NumOfSeries)[1:]):
            fit_params['zpfit_%i' % listOfbootsnoFirst[iy]].expr='zpfit_%i'% listOfboots[FirstSingularitem]

        for iy in tuple(range(len(listOfbootsnoFirst))):#tuple(range(NumOfSeries)[1:]):
            fit_params['elevationfit_%i' % listOfbootsnoFirst[iy]].expr='elevationfit_%i'% listOfboots[FirstSingularitem]

        for iy in tuple(range(len(listOfbootsnoFirst))):#tuple(range(NumOfSeries)[1:]):
            fit_params['cfit_%i' % listOfbootsnoFirst[iy]].expr='cfit_%i'% listOfboots[FirstSingularitem]

        mini = Minimizer(objective_diffraction, fit_params,nan_policy = "raise",scale_covar=False,iter_cb=per_iteration, fcn_args=(zetaprime2,dataBooted,1.),)

        result = mini.leastsq(ftol=1.e-2)

        report_fit(result.params,show_correl=True, min_correl=0.01)

        #residLst.append(result.residual)
        residLst = result.residual
        
        residsqrlst.append(np.sum(residLst**2))

        UnDi.SingleAmplitudeDone = 0

        #exec(open("./saving_things.py").read())

        #if len(WLidxList) > 1:
        #    deltaGuess = result.params['delta0_%i'% listOfboots[FirstSingularitem]].value
        #    print("Inherit delta0fit")



        exec(open("./plot_modular_results_NonGUI.py").read())


#ModFit(zetaprime2,List825VolumeMinusBG4Mean,WidthOfStripe,pixelof404nm,dispersion,CalibWl,NumOfSeries,LengthOfData,Interval,[263],list(range(200))) ##Position variants
"""
def ModFit(
zetaprime2,
List825VolumeMinusBG4Mean,
WidthOfStripe,
pixelof404nm,
dispersion,
CalibWl,
NumOfSeries,
LengthOfData,
Interval,
WLidxList,
PositionidxList):
"""
ModFit(zetaprime2,List825VolumeMinusBG4Mean,WidthOfStripe,pixelof404nm,dispersion,CalibWl,NumOfSeries,LengthOfData,Interval,list(range(10,460,10)),[0])#,list(range(1,460,1)),[0]) ##WL variants

exit()