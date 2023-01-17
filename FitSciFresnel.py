import matplotlib.pyplot as plt

import UndulatorDiffractionOneDimensional as UnDi

import AmplitudeFormula as AUndu

from numpy import transpose as nptp
import numpy as np

import cv2
import os , sys
import re

import csv

from natsort import natsorted
from scipy.optimize import curve_fit
from lmfit import minimize, Parameters, report_fit, Minimizer

ratio = 1000. / 606.#606. #how many pixel are aperture
UnDi.apy = .004


zetaprime2 = np.linspace(-ratio*UnDi.apy,ratio*UnDi.apy,1000)
fixedsize = .245
k0 = 0
k = k0

NumOfStripes = 1
NumOfSeries = 825

def per_iteration(pars, iteration, resid, *args, **kws):
    print(" resid ", sum(resid))#sum(abs(np.array(resid)))/1000.)
    #print(" resid ", resid)
    print("Amplfit: " + str(pars['Amplfit_1'].value))
    print("bfit: " + str(pars['bfit_1'].value))
    print('Ytran: ' + str(pars['Ytran_1'].value))
    print("gammafit: " + str(pars['gamma_1'].value))
    print("delta0fit: " + str(pars['delta0_1'].value))

    Amplfit = pars['Amplfit_1'].value
    Ytran = pars['Ytran_1'].value
    gammafit = pars['gamma_1'].value
    delta0fit = pars['delta0_1'].value
    bfit = pars['bfit_1'].value
    residLst = []
    
    UnDi.SingleAmplitudeDone = 0
    UnDi.AmplXY = []
    #cv2.imshow("resid",nptp(abs(resid).reshape((200, 1000))))
    #cv2.waitKey(0)
    
    #if sum(resid)>-7214.:
    #    return True
    
    #dataLst = []
    #print(len(datasigma))
    #print(datasigma[0])
    
    #for Strpidx in range(NumOfSeries):
    #    residLst.append(np.array((UnDi.diffractionintensitySci(zetaprime2,Amplfit,Ytran,delta0fit+.001*Strpidx,gammafit)-data[Strpidx])/datasigma[Strpidx]))
    #print(" My resid ", sum(np.array(residLst).flatten()))#sum(abs(np.array(residLst)))/1000.)
    #exit()

    

def diffraction_dataset(params,i,yy):
    Amplfit = params['Amplfit_%i' % (i+1)].value
    Ytran   = params['Ytran_%i' % (i+1)].value
    gammafit  = params['gamma_%i' % (i+1)].value
    delta0fit  = params['delta0_%i' % (i+1)].value
    bfit  = params['bfit_%i' % (i+1)].value
    #print(delta0fit+i*.001)
    if i > 0:
        UnDi.SingleAmplitudeDone = 1
    Diff_intensity = UnDi.diffractionintensity(yy,Ytran,i,AUndu.z0,404.3e-9,delta0fit+i*.001,gammafit,Amplfit,bfit)


    return Diff_intensity

def objective_diffraction(params, yy, data,weights):
    ndata,nx = data.shape
    
    resid = 0.0*data[:]
    
    for Strpidx in range(ndata):
        resid[Strpidx] = (diffraction_dataset(params,Strpidx,yy) - data[Strpidx])*weights[Strpidx]#**2.
        
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

slashes = [i for i,ltr in enumerate(DataPath) if ltr == "\\"]
DataPathBase = DataPath[0:slashes[-2]]

print(DataPathBase)

#datarangelow = 1804 + 550

Backgroundimage = cv2.imread(DataPathBase + "/BG/raw_shortExp/File1.tif",cv2.IMREAD_UNCHANGED)
ListOfAllMatrices = []
ListOfAllSigmaMatrices = []


WidthOfStripe = 5#int((2304-datarangelow)/NumOfStripes)
print(WidthOfStripe)

MatrixOfSlit = []

maybeBadGamma = 0

targetfiles = ["Amplfit.csv",
               "Ytran.csv",
               "gammafit.csv",
               "delta0fit.csv",
               "bfit.csv",
               "Amplfitstderr.csv",
               "Ytranstderr.csv",
               "gammafitstderr.csv",
               "delta0fitstderr.csv",
               "bfitstderr.csv",
               "amplYtranCorrel.csv",
               "amplgammafitCorrel.csv",
               "ampldelta0fitCorrel.csv",
               "amplbfitCorrel.csv",
               "YtrangammafitCorrel.csv",
               "Ytrandelta0fitCorrel.csv",
               "YtranbfitCorrel.csv",
               "gammadelta0fitCorrel.csv",
               "gammabfitCorrel.csv",
               "delta0fitbfitCorrel.csv"]
badk = []

DoNotSave = 0

while k < int(len(filelist)/4):
    #Guessvalues = [4.4,-.0004,delta0+k*.001]
    if k == k0:
        checklist = os.listdir(DataPath + "\SharedP")
        maybeBadGamma = 0
        for ckid in range(len(targetfiles)):
            if targetfiles[ckid] in checklist:
                os.remove(DataPath + "\SharedP" + "//" + targetfiles[ckid])

    SigmaOneLineEachFromTenPicturesAtOnePosition=[]

    SeriesOneLineEachFromFourPicturesAtOnePosition=[]
    SeriesOneLineEachFromFourPicturesAtOnePositionForSigma=[]

    for ii in range(NumOfSeries):
        OneLineEachFromFourPicturesAtOnePosition = []
        OneLineEachFromFourPicturesAtOnePositionForSigma=[]
        for i in range(4):
            noti = i
            PictureMinusBG = cv2.imread(DataPath + "/" + filelist[i+ii*4+k*4],cv2.IMREAD_UNCHANGED).astype(int)-Backgroundimage
            PictureMinusBGForSigma = cv2.imread(DataPath + "/" + filelist[i+ii*4+k*4],cv2.IMREAD_UNCHANGED).astype(int)-Backgroundimage

            PictureMinusBG = np.float32(PictureMinusBG)
            PictureMinusBGForSigma = np.float32(PictureMinusBGForSigma)
            #print(len(PictureMinusBG[0]))

            PictureMinusBG = nptp(nptp(PictureMinusBG)[(2304-WidthOfStripe):])
            #print(len(PictureMinusBG[0]))
            PictureMinusBGForSigma = nptp(nptp(PictureMinusBGForSigma)[(2304-WidthOfStripe):])

            OneLineEachFromFourPicturesAtOnePosition.append(PictureMinusBG)
            OneLineEachFromFourPicturesAtOnePositionForSigma.append(PictureMinusBGForSigma)
            
        
        #print(len(OneLineEachFromFourPicturesAtOnePositionForSigma))
        #print(len(OneLineEachFromFourPicturesAtOnePositionForSigma[0]))
        #print(len(OneLineEachFromFourPicturesAtOnePositionForSigma[0][0]))
        
        TenPixelFromTenPicturesTimesLineAtOnePositionForSigma = np.std(OneLineEachFromFourPicturesAtOnePositionForSigma, axis=0)
        
        #print(len(OneLineEachFromFourPicturesAtOnePositionForSigma))
        #print(len(OneLineEachFromFourPicturesAtOnePositionForSigma[0]))
        
        #cv2.imshow("imahe",OneLineEachFromFourPicturesAtOnePositionForSigma)
        #cv2.waitKey(0)
        
        TenPixelFromTenPicturesTimesLineAtOnePosition = sum(OneLineEachFromFourPicturesAtOnePosition)/4.
        #TenPixelFromTenPicturesTimesLineAtOnePositionForSigma = sum(OneLineEachFromFourPicturesAtOnePositionForSigma)/4.
        
        SeriesOneLineEachFromFourPicturesAtOnePosition.append(TenPixelFromTenPicturesTimesLineAtOnePosition)
        SeriesOneLineEachFromFourPicturesAtOnePositionForSigma.append(TenPixelFromTenPicturesTimesLineAtOnePositionForSigma)

    print(len(SeriesOneLineEachFromFourPicturesAtOnePosition))
    print(len(SeriesOneLineEachFromFourPicturesAtOnePosition[0]))
    print(len(SeriesOneLineEachFromFourPicturesAtOnePosition[0][0]))

    interpolatedpicture = []
    interpolatedsigmapicture = []
    
    for ii in range(NumOfSeries):
        stripe = SeriesOneLineEachFromFourPicturesAtOnePosition[ii]#nptp(nptp(SeriesOneLineEachFromFourPicturesAtOnePosition[ii])[WidthOfStripe*0:WidthOfStripe+WidthOfStripe*0])
        stripesigma = SeriesOneLineEachFromFourPicturesAtOnePositionForSigma[ii]#nptp(nptp(SeriesOneLineEachFromFourPicturesAtOnePositionForSigma[ii])[WidthOfStripe*0:WidthOfStripe+WidthOfStripe*0])
        
        #cv2.imshow("imahe",stripesigma)
        #cv2.waitKey(0)

        stripe = cv2.resize(stripe,(1,len(TenPixelFromTenPicturesTimesLineAtOnePosition)),interpolation=cv2.INTER_AREA) #1 zu 1000
        stripesigma = cv2.resize(stripesigma,(1,len(TenPixelFromTenPicturesTimesLineAtOnePositionForSigma)),interpolation=cv2.INTER_AREA) #1 zu 1000

        interpolatedpicture.append(nptp(stripe)[0])
        interpolatedsigmapicture.append(nptp(stripesigma))
        
        #print(interpolatedpicture)
        #print(interpolatedsigmapicture)
        
        #data2D = np.transpose(stripe)[0]/200.
        #datasigma2D = np.transpose(stripesigma)[0]/200.
        
        print(stripe[500])

    data = np.array(interpolatedpicture)
    datasigma = np.array(interpolatedsigmapicture)


    dataLst = []
    for Strpidx in range(NumOfSeries):
        dataLst.append(data[Strpidx][500])
    """
    fntsze = 13
    figA = plt.figure("Diffraction")
    axA = figA.add_subplot(111)
    axA.plot(range(len(dataLst)), dataLst, '-')
    axA.set_xlabel('Position on Aperture [mm]',fontsize=fntsze)
    axA.set_ylabel('Intensity [a.u.]',fontsize=fntsze)
    axA.tick_params(axis='both', which='major', labelsize=fntsze)

    plt.show()
    """


    print(k)
    print((NumOfSeries/5.))
    
    Guessvalueslmfit = [3.2,.0001,351.0,.195+k*.001,1.]
    #Guessvalueslmfit = [5.28,-.00037,351.941,.196773+k*.001,2.05]

    fit_params = Parameters()
    for iy,y, in enumerate(data):
        fit_params.add('Amplfit_%i' % (iy+1),value=Guessvalueslmfit[0])#, min = Guessvalueslmfit[0]-1.1, max=Guessvalueslmfit[0]+1.1)
        fit_params.add('Ytran_%i' % (iy+1),value=Guessvalueslmfit[1])#, min = Guessvalueslmfit[1]-.0001, max=Guessvalueslmfit[1]+.0001)
        fit_params.add('gamma_%i' % (iy+1),value=Guessvalueslmfit[2])#,min = Guessvalueslmfit[2]-1., max = Guessvalueslmfit[2]+1.)
        fit_params.add('delta0_%i' % (iy+1),value=Guessvalueslmfit[3])#, min = Guessvalueslmfit[3]-.003, max=Guessvalueslmfit[3]+.003)
        fit_params.add('bfit_%i' % (iy+1),value=Guessvalueslmfit[4])#, min = Guessvalueslmfit[4]-.1, max=Guessvalueslmfit[4]+.1)

    for iy in tuple(range(NumOfSeries+1)[2:]):
        fit_params['Amplfit_%i' % iy].expr='Amplfit_1'
        
    for iy in tuple(range(NumOfSeries+1)[2:]):
        fit_params['Ytran_%i' % iy].expr='Ytran_1'
        
    for iy in tuple(range(NumOfSeries+1)[2:]):
        fit_params['gamma_%i' % iy].expr='gamma_1'
        
    for iy in tuple(range(NumOfSeries+1)[2:]):
        fit_params['delta0_%i' % iy].expr='delta0_1'

    for iy in tuple(range(NumOfSeries+1)[2:]):
        fit_params['bfit_%i' % iy].expr='bfit_1'

    mini = Minimizer(objective_diffraction, fit_params,nan_policy = "raise",scale_covar=False,iter_cb=per_iteration, fcn_args=(zetaprime2,data,1./datasigma),)

    #print(data[0][500])
    #print(data[1][500])

    print(len(zetaprime2))
    print(len(datasigma[0]))


    result = mini.leastsq(ftol=1.e-2)
    """
    report_fit(result.params,show_correl=True, min_correl=0.01)

    y_fitLst = []
    dataLst = []
    


    for Strpidx in range(NumOfSeries):

        WLofPixel = 404.184e-9# -.478e-9*1.#*(Strpidx-(int(NumOfStripes/2)))
        Sourcelst=[]

        Amplfit = result.params['Amplfit_1'].value
        Ytran = result.params['Ytran_1'].value
        gammafit = result.params['gamma_1'].value
        delta0fit = result.params['delta0_1'].value
        bfit = result.params['bfit_1'].value

        Sourcelst = UnDi.plainintensity(zetaprime2,Ytran,Strpidx,AUndu.z0,404.3e-9,delta0fit+.001*Strpidx,gammafit,Amplfit,bfit)
        
        #y_fitLst.append(sum(UnDi.diffractionintensity(zetaprime2,Ytran,Strpidx,AUndu.z0,404.3e-9,delta0fit+Strpidx*.001,gammafit,Amplfit,bfit)[475:485]))
        #dataLst.append(sum((data[Strpidx])[475:485]))
        y_fitLst.append(UnDi.diffractionintensity(zetaprime2,Ytran,Strpidx,AUndu.z0,404.3e-9,delta0fit+Strpidx*.001,gammafit,Amplfit,bfit)[480])
        dataLst.append((data[Strpidx])[480])

    fntsze = 13
    figA = plt.figure("Diffraction")
    axA = figA.add_subplot(111)

    #axA.annotate(r'$\Delta d + d = %d $' % k, xy=(0, maxofplt), fontsize=fntsze)
    axA.plot(range(len(dataLst)), dataLst, '-', range(len(dataLst)), y_fitLst, '-')
    #axA.set_xlim(-7, 7)
    #axA.set_ylim(-10, maxofplt+10)
    axA.set_xlabel('Position on Aperture [mm]',fontsize=fntsze)
    axA.set_ylabel('Intensity [a.u.]',fontsize=fntsze)
    axA.tick_params(axis='both', which='major', labelsize=fntsze)
    plt.show()
    
    axis_font = {'fontname':'Arial', 'size':'21.5'}
    numbers_font_Century = {'fontname':'Arial', 'size':'21.5'}
    
    y_fitLst = []

    for Strpidx in range(NumOfSeries):

        #y_fit = diffraction_dataset(result.params, Strpidx, zetaprime2)
        #y_fit = UnDi.diffractionintensitySci(zetaprime2,*popt)#
        #y_fit = UnDi.diffractionintensitySci_delta(zetaprime2,deltad)#diffraction_dataset(result.params, Strpidx, zetaprime2)

        WLofPixel = 404.184e-9# -.478e-9*1.#*(Strpidx-(int(NumOfStripes/2)))
        Sourcelst=[]
        
        Amplfit = result.params['Amplfit_1'].value
        Ytran = result.params['Ytran_1'].value
        gammafit = result.params['gamma_1'].value
        delta0fit = result.params['delta0_1'].value
        bfit = result.params['bfit_1'].value

        #deltad = delta0+Strpidx*.001

        Sourcelst = UnDi.plainintensity(zetaprime2,Ytran,Strpidx,AUndu.z0,404.3e-9,delta0fit+.001*Strpidx,gammafit,Amplfit,bfit)

        SourcelstIntensity1D = 55000.*abs(np.array(Sourcelst))

        y_fit = UnDi.diffractionintensity(zetaprime2,Ytran,Strpidx,AUndu.z0,404.3e-9,delta0fit+Strpidx*.001,gammafit,Amplfit,bfit)
               #UnDi.diffractionintensitySci(zetaprime2,Ytran,delta0fit+.001*Strpidx,gammafit,Amplfit,bfit)

        #residualbeforeplot = sum(((np.array(y_fit-data[Strpidx, :])[0:1000])**2.)/(1000.*datasigma[Strpidx][0:1000]**2.))
        #print(residualbeforeplot)

        #innerresidual = sum(((np.array(y_fit-data[Strpidx, :])[400:600])**2.)/(200.*datasigma[Strpidx][400:600]**2.))
        #print(innerresidual)

        #maxofdata = np.around(max(data[Strpidx, :]), decimals=-1)
        #maxofSourcelstIntensity1D = np.around(max(SourcelstIntensity1D), decimals=-1)
        #maxofplt = max([maxofdata,maxofSourcelstIntensity1D])
        
        y_fitLst.append(y_fit)

        fntsze = 13

        figA = plt.figure("Diffraction")
        axA = figA.add_subplot(111)
        figB = plt.figure("Sources")
        axB = figB.add_subplot(111)

        #axA.annotate(r'$\Delta d + d = %d $' % k, xy=(0, maxofplt), fontsize=fntsze)
        axA.plot(1000.*zetaprime2[0:1000], data[Strpidx, :], '-', 1000.*zetaprime2, y_fit, '-')
        axA.set_xlim(-7, 7)
        #axA.set_ylim(-10, maxofplt+10)
        axA.set_xlabel('Position on Aperture [mm]',fontsize=fntsze)
        axA.set_ylabel('Intensity [a.u.]',fontsize=fntsze)
        axA.tick_params(axis='both', which='major', labelsize=fntsze)
        
        #axB.annotate(r'$\Delta d + d = %d $' % k, xy=(0, maxofplt), fontsize=fntsze)
        axB.plot(1000.*zetaprime2, SourcelstIntensity1D, '-')
        axB.set_xlim(-7, 7)
        #axB.set_ylim(-10, maxofplt+10)
        axB.set_xlabel('Position on Aperture [mm]',fontsize=fntsze)
        axB.set_ylabel('Intensity [a.u.]',fontsize=fntsze)
        axB.tick_params(axis='both', which='major', labelsize=fntsze)

        plt.show()

    #"""
    
    Amplfit = result.params['Amplfit_1'].value
    Ytran = result.params['Ytran_1'].value
    gammafit = result.params['gamma_1'].value
    delta0fit = result.params['delta0_1'].value
    bfit = result.params['bfit_1'].value

    Amplfitstderr = result.params['Amplfit_1'].stderr
    Ytranstderr = result.params['Ytran_1'].stderr
    gammafitstderr = result.params['gamma_1'].stderr
    delta0fitstderr = result.params['delta0_1'].stderr
    bfitstderr = result.params['bfit_1'].stderr

    amplYtranCorrel = result.covar[0][1] / (result.params["Amplfit_1"].stderr*result.params["Ytran_1"].stderr)
    amplgammafitCorrel = result.covar[0][2] / (result.params["Amplfit_1"].stderr*result.params["gamma_1"].stderr)
    ampldelta0fitCorrel = result.covar[0][3] / (result.params["Amplfit_1"].stderr*result.params["delta0_1"].stderr)
    amplbfitCorrel = result.covar[0][4] / (result.params["Amplfit_1"].stderr*result.params["bfit_1"].stderr)
    YtrangammafitCorrel = result.covar[1][2] / (result.params["Ytran_1"].stderr*result.params["gamma_1"].stderr)
    Ytrandelta0fitCorrel = result.covar[1][3] / (result.params["Ytran_1"].stderr*result.params["delta0_1"].stderr)
    YtranbfitCorrel = result.covar[1][4] / (result.params["Ytran_1"].stderr*result.params["bfit_1"].stderr)
    gammadelta0fitCorrel = result.covar[2][3] / (result.params["gamma_1"].stderr*result.params["delta0_1"].stderr)
    gammabfitCorrel = result.covar[2][4] / (result.params["gamma_1"].stderr*result.params["bfit_1"].stderr)
    delta0fitbfitCorrel = result.covar[3][4] / (result.params["delta0_1"].stderr*result.params["bfit_1"].stderr)    

    if DoNotSave == 0:
        print("saving")
        with open(DataPath + "\\SharedP\\Amplfit.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([Amplfit])

        with open(DataPath + "\\SharedP\\Ytran.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([Ytran])

        with open(DataPath + "\\SharedP\\gammafit.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([gammafit])

        with open(DataPath + "\\SharedP\\delta0fit.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([delta0fit])

        with open(DataPath + "\\SharedP\\bfit.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([bfit])

        with open(DataPath + "\\SharedP\\Amplfitstderr.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([Amplfitstderr])

        with open(DataPath + "\\SharedP\\Ytranstderr.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([Ytranstderr])

        with open(DataPath + "\\SharedP\\gammafitstderr.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([gammafitstderr])

        with open(DataPath + "\\SharedP\\delta0fitstderr.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([delta0fitstderr])

        with open(DataPath + "\\SharedP\\bfitstderr.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([bfitstderr])



        with open(DataPath + "\\SharedP\\amplYtranCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([amplYtranCorrel])

        with open(DataPath + "\\SharedP\\amplgammafitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([amplgammafitCorrel])
        
        with open(DataPath + "\\SharedP\\ampldelta0fitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([ampldelta0fitCorrel])
            
        with open(DataPath + "\\SharedP\\amplbfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([amplbfitCorrel])

        with open(DataPath + "\\SharedP\\YtrangammafitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([YtrangammafitCorrel])

        with open(DataPath + "\\SharedP\\Ytrandelta0fitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([Ytrandelta0fitCorrel])

        with open(DataPath + "\\SharedP\\YtranbfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([YtranbfitCorrel])
        
        with open(DataPath + "\\SharedP\\gammadelta0fitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([gammadelta0fitCorrel])
            
        with open(DataPath + "\\SharedP\\gammabfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([gammabfitCorrel])

        with open(DataPath + "\\SharedP\\delta0fitbfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([delta0fitbfitCorrel])
    exit()