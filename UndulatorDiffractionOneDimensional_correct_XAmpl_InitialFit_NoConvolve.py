import matplotlib.pyplot as plt

import InterferenceFunction_more as IFunc
import AmplitudeFormula as AUndu

from scipy import ndimage

from numpy import transpose as nptp
import numpy as np

global apy
global SingleAmplitudeDone
global AmplXY
global CMOS

#apy = .0065/2.*.16

SingleAmplitudeDone = 0
gfdsgfdsg = 1
AmplXY = []

#apy = .014
#apx = .004

def ConditionAmpl(yy,xx,Yscrfit,Yt):
    if xx == 0:
        if (yy-Yscrfit-Yt) < 0.:
            phi = 2.*np.pi-np.pi/2.
        else:
            phi = np.pi/2.
    else:
        if yy>0:
            if xx>0:
                phi = np.arctan(yy/xx)
            else:
                phi = np.pi + np.arctan(yy/xx)
        else:
            if xx>0:
                phi = 2.*np.pi + np.arctan(yy/xx)
            else:
                phi = np.pi + np.arctan(yy/xx)
    return phi


def Conditionphase(yy,xx,Yscrfit,Yt):
    if xx == 0:
        if (yy-Yscrfit-Yt) < 0.:
            phi = 2.*np.pi-np.pi/2.
        else:
            phi = np.pi/2.
    else:
        if yy>0:
            if xx>0:
                phi = np.arctan(yy/xx)
            else:
                phi = np.pi + np.arctan(yy/xx)
        else:
            if xx>0:
                phi = 2.*np.pi + np.arctan(yy/xx)
            else:
                phi = np.pi + np.arctan(yy/xx)
    return phi


def amplitudeonaperture_X(yy,xx,Ytran,WLofPixel,LambdaRfit,z0,deltad,gammafit,Yscreenfit): # Xtran=Xtranslation
    Sourcelst = []
    Ku = np.sqrt(abs(2.*gammafit**2*LambdaRfit/AUndu.LambU-1.)*2.)

    Sourcelst = [[AUndu.Amplitude(gammafit,WLofPixel,LambdaRfit,Ku,(xx[i]**2. + (yy[j]-Yscreenfit-Ytran)**2.)**.5/z0,ConditionAmpl(yy[0],xx[i],Yscreenfit,Ytran)) for j in range(len(yy))] for i in range(len(xx))]
    return Sourcelst

def phaseonaperture_X(yy,xx,Ytran,LambdaRfit,z0,deltad,gammafit,Yscreenfit):
    phaselst = []
    Ku = np.sqrt(abs(2.*gammafit**2*LambdaRfit/AUndu.LambU-1.)*2.)

    phaselst = [[deltad/(2.*gammafit**2)*(1+((xx[i]**2. + (yy[j]-Yscreenfit-Ytran)**2.)**.5/z0)**2*gammafit**2) + AUndu.Lu/(4.*gammafit**2)*(2.+Ku**2 +2.*gammafit**2*((xx[i]**2. + (yy[j]-Yscreenfit-Ytran)**2.)**.5/z0)**2) for j in range(len(yy))] for i in range(len(xx))]
    return phaselst

def amplitudeonaperture(yy,Ytran,WLofPixel,LambdaRfit,z0,deltad,gammafit,Yscreenfit): # Xtran=Xtranslation
    Sourcelst = []
    Ku = np.sqrt(abs(2.*gammafit**2*LambdaRfit/AUndu.LambU-1.)*2.)

    for i in range(len(yy)):
        thetaradius = abs(yy[i]-Yscreenfit-Ytran) #20200319 -5.e-4 18px
        theta = thetaradius/z0

        if (yy[i]-Yscreenfit-Ytran) < 0.:
            phi = -np.pi/2.
        else:
            phi = np.pi/2.
        #if i ==1:
        #    print(Ampl)

        #print(AUndu.Amplitude(gammafit,WLofPixel,LambdaRfit,Ku,theta,phi)*np.array([np.exp(1.j * 2.*np.pi/WLofPixel),np.exp(1.j * 2.*np.pi/WLofPixel)]))
        #exit()
        Ampl = AUndu.Amplitude(gammafit,WLofPixel,LambdaRfit,Ku,theta,phi)#*np.array([np.sqrt(1.+deltad/10),np.sqrt(1.+deltad/10)])
        Sourcelst.append(Ampl)
    return Sourcelst

def phaseonaperture(yy,Ytran,LambdaRfit,z0,deltad,gammafit,Yscreenfit):
    phaselst = []
    Ku = np.sqrt(abs(2.*gammafit**2*LambdaRfit/AUndu.LambU-1.)*2.)

    for i in range(len(yy)):
        thetaradius = abs(yy[i]-Yscreenfit-Ytran) #20200319 -5.e-4 18px
        theta = thetaradius/z0

        if (yy[i]-Yscreenfit-Ytran) < 0.:
            phi = -np.pi/2.
        else:
            phi = np.pi/2.
        
        phaseR =  deltad/(2.*gammafit**2)*(1+theta**2*gammafit**2) + AUndu.Lu/(4.*gammafit**2)*(2.+Ku**2 +2.*gammafit**2*theta**2)
        phaselst.append(phaseR)
    return phaselst

def plotaperture(zetaprime2,Sourcelst2D):
    fig2 = plt.figure("Original field nr1")
    plt.xlabel('horizontal')
    plt.ylabel('vertical')
    plt.plot(1000.*zetaprime2[0:int(len(zetaprime2)/2.)*2],Sourcelst2D)
    plt.show()

def superamplitudeInitialFit(yy,Ytran,Strpidx,z0,WLofPixel,LambdaRfit,deltad,gammafit,Amplfit,bfit,Yscreenfit):

    if SingleAmplitudeDone == 0:
        AmplXY.append(np.array(nptp(amplitudeonaperture(yy,Ytran,WLofPixel,LambdaRfit,z0,deltad,gammafit,Yscreenfit))))
        print("idx" ,Strpidx)

    phaseR = phaseonaperture(yy,Ytran,LambdaRfit,z0,deltad,gammafit,Yscreenfit)

    AmplTwoX = []
    AmplTwoY = []
    bfitTwo  = []

    for i in range(len(yy)):
        if yy[i] > -apy+Yscreenfit:
            if yy[i] < apy+Yscreenfit:
                bfitTwo.append(bfit)
            else:
                bfitTwo.append(.0)
                #bfitTwo.append(.7)
                #bfitTwo.append(.4)
        else:
            bfitTwo.append(.0)

        if (yy[i] > -apy+Yscreenfit and yy[i] < apy+Yscreenfit):
            #print(yy[i]-Ytran)
            AmplX,AmplY = AmplXY[0][0],AmplXY[0][1]
            AmplTwoX.append(AmplX[i]*Amplfit)#*(1.+np.exp(1.j*2.*np.pi/WLofPixel*phaseR[i]))+bfit)
            AmplTwoY.append(AmplY[i]*Amplfit)#*(1.+np.exp(1.j*2.*np.pi/WLofPixel*phaseR[i]))+bfit)

        else:
            AmplTwoX.append(0.)
            AmplTwoY.append(0.)
    return [np.array(AmplTwoX),np.array(AmplTwoY),np.array(phaseR),np.array(bfitTwo)]

def superamplitudeInitialFit_X(yy,xx,Ytran,Strpidx,z0,WLofPixel,LambdaRfit,deltad,gammafit,Amplfit,bfit,Yscreenfit):
    #AmplXY = []
    if SingleAmplitudeDone == 0:
        AmplXY.append(np.array(nptp(amplitudeonaperture_X(yy,xx,Ytran,WLofPixel,LambdaRfit,z0,deltad,gammafit,Yscreenfit))))
        print("idx" ,Strpidx)
    #print("AmplXY",AmplXY)
    #phaseonaperture(yy,Ytran,LambdaRfit,z0,deltad,gammafit,Yscreenfit)

    phaseR = phaseonaperture_X(yy,xx,Ytran,LambdaRfit,z0,deltad,gammafit,Yscreenfit)

    #print("phaseR",np.array(phaseR))

    #print(AmplXY)
    #exit()
    #AmplXY = [np.transpose(AmplXY)]
    AmplTwoX = []
    AmplTwoY = []
    bfitTwo  = []
    for i in range(len(xx)):
        bfitTwo1D = []
        AmplTwoX1D = []
        AmplTwoY1D = []
        for j in range(len(yy)):

            if yy[j] > -apy+Yscreenfit:
                if yy[j] < apy+Yscreenfit:
                    bfitTwo1D.append(bfit)
                else:
                    bfitTwo1D.append(.0)
                    #bfitTwo.append(.7)
                    #bfitTwo.append(.4)
            else:
                bfitTwo1D.append(.0)

            if (yy[j] > -apy+Yscreenfit and yy[j] < apy+Yscreenfit):
                AmplX,AmplY = AmplXY[0][0],AmplXY[0][1]

                AmplTwoX1D.append(AmplX[j][i]*Amplfit)#*(1.+np.exp(1.j*2.*np.pi/WLofPixel*phaseR[i]))+bfit)
                AmplTwoY1D.append(AmplY[j][i]*Amplfit)#*(1.+np.exp(1.j*2.*np.pi/WLofPixel*phaseR[i]))+bfit)

            else:
                AmplTwoX1D.append(np.array(0.))
                AmplTwoY1D.append(np.array(0.))
                
                

        AmplTwoX.append(AmplTwoX1D)
        AmplTwoY.append(AmplTwoY1D)
        bfitTwo.append(bfitTwo1D)

    bfitTwo = np.array(bfitTwo)
    AmplTwoX = np.array(AmplTwoX)
    AmplTwoY = np.array(AmplTwoY)
    phaseR = np.array(phaseR)

    return [np.array(AmplTwoX),np.array(AmplTwoY),np.array(phaseR),np.array(bfitTwo)]

def plainintensity(yy,xx,Ytran,Strpidx,z0,WLofPixel,LambdaRfit,deltad,gammafit,Amplfit,bfit,zpfit,Yscreenfit,elevationfit):

    AmplTwoU1 = superamplitudeInitialFit_X(yy,xx,Ytran,Strpidx,z0,WLofPixel,LambdaRfit,deltad,gammafit,Amplfit,bfit,Yscreenfit)
    AmplTwoXU1, AmplTwoYU1 = AmplTwoU1[0], AmplTwoU1[1]

    AmplTwoU2 = superamplitudeInitialFit_X(yy,xx,Ytran,Strpidx,z0-deltad,WLofPixel,LambdaRfit,deltad,gammafit,Amplfit,bfit,Yscreenfit)

    phaseR = AmplTwoU2[2]
    bfitTwo = AmplTwoU2[3]

    #phaseR = np.transpose(phaseR)
    
    #AmplTwoXU2, AmplTwoYU2 = AmplTwoU2[0]*np.exp(1.j*2.*np.pi/WLofPixel*phaseR) + bfitTwo, AmplTwoU2[1]*np.exp(1.j*2.*np.pi/WLofPixel*phaseR) + bfitTwo
    AmplTwoXU2, AmplTwoYU2 = AmplTwoU2[0]*np.exp(1.j*2.*np.pi/WLofPixel*phaseR) + bfitTwo, AmplTwoU2[1]*np.exp(1.j*2.*np.pi/WLofPixel*phaseR) + bfitTwo

    #ElectronSource = np.array(IFunc.ElectronSourcefunction(yy))

    #AmplTwoXU1 = IFunc.convolveplain(AmplTwoXU1,ElectronSource)
    #AmplTwoXU2 = IFunc.convolveplain(AmplTwoXU2,ElectronSource)
    #AmplTwoYU1 = IFunc.convolveplain(AmplTwoYU1,ElectronSource)
    #AmplTwoYU2 = IFunc.convolveplain(AmplTwoYU2,ElectronSource)

    Sourcelst1D1000 = np.conjugate(AmplTwoXU1+AmplTwoXU2)*(AmplTwoXU1+AmplTwoXU2)#+np.conjugate(AmplTwoYU1+AmplTwoYU2)*(AmplTwoYU1+AmplTwoYU2)

    Sourcelst1D1000 = np.real(Sourcelst1D1000)

    #Sourcelst1D1000 = sum(np.sqrt(Sourcelst1D1000))

    return Sourcelst1D1000*10e7*(1+elevationfit*deltad/100)

def diffractionintensityInitialFit_X(yy,xx,Ytran,Strpidx,z0,WLofPixel,LambdaRfit,deltad,gammafit,Amplfit,bfit,zpfit,Yscreenfit,elevationfit,whichdata,cfit):

    AmplTwoU1 = superamplitudeInitialFit_X(yy,xx,Ytran,Strpidx,z0,WLofPixel,LambdaRfit,deltad,gammafit,Amplfit,bfit,Yscreenfit)
    AmplTwoXU1, AmplTwoYU1 = AmplTwoU1[0][0], AmplTwoU1[1][0]

    AmplTwoU2 = superamplitudeInitialFit_X(yy,xx,Ytran,Strpidx,z0-deltad,WLofPixel,LambdaRfit,deltad,gammafit,Amplfit,bfit,Yscreenfit)
    phaseR = AmplTwoU2[2][0]
    bfitTwo = AmplTwoU2[3][0]

    AmplTwoXU2, AmplTwoYU2 = AmplTwoU2[0][0], AmplTwoU2[1][0]

    Sourcelst1D1000_1 = AmplTwoXU1
    Sourcelst1D1000_2 = AmplTwoXU2

    ElectronSource = np.array(IFunc.ElectronSourcefunction(yy))

    fresnelforward1 = np.array(IFunc.squaredphase(yy,Ytran*0.,0,zpfit,WLofPixel)) #Here zero is true, it is seen as complementary to deltad + Lu
    fresnelforward2 = np.array(IFunc.squaredphase(yy,Ytran*0.,deltad + AUndu.Lu,zpfit,WLofPixel))

    origprime1 = IFunc.convolvefresnel(Sourcelst1D1000_1,fresnelforward1,ElectronSource)/(WLofPixel*100.)
    origprime2 = IFunc.convolvefresnel(Sourcelst1D1000_2*np.exp(1.j*2.*np.pi/WLofPixel*phaseR),fresnelforward2,ElectronSource)/(WLofPixel*100.) + bfitTwo

    if whichdata == 1:
        Sourcelst1D1000 = (np.conjugate(origprime1+origprime2)*(origprime1+origprime2))/20.
        Sourcelst1D1000 = Sourcelst1D1000.real
        Sourcelst1D1000 = Sourcelst1D1000*(1+elevationfit*deltad/100) - np.array([.05]*240+[cfit]*606+[.005]*(1000-606-240)) # This is normal
        #Sourcelst1D1000 = Sourcelst1D1000*(1+elevationfit*deltad/100) - np.array([.05]*232+[cfit]*606+[.005]*(1000-606-232))
        #print("Is NoConvolveData_825.py")
    else:
        Sourcelst1D1000 = (np.conjugate(origprime1+origprime2)*(origprime1+origprime2))/20.
        Sourcelst1D1000 = Sourcelst1D1000.real
        Sourcelst1D1000 = Sourcelst1D1000*(1+elevationfit*deltad/100) - np.array([cfit]*240+[cfit]*606+[cfit]*(1000-606-240))
        #print("Is NoConvolveData_825_SplitYTran.py")


    #Sourcelst1D1000 = Sourcelst1D1000*(1+elevationfit*deltad/100) - np.array([1.]*240+[cfit]*606+[.5]*(1000-606-240))
    #Sourcelst1D1000 = Sourcelst1D1000*(1+elevationfit*deltad/100) - np.array([.05]*240+[cfit]*606+[.005]*(1000-606-240)) # This is normal


    return Sourcelst1D1000
    
def diffractionintensityInitialFit(yy,xx,Ytran,Strpidx,z0,WLofPixel,LambdaRfit,deltad,gammafit,Amplfit,bfit,zpfit,Yscreenfit,elevationfit):

    AmplTwoU1 = superamplitudeInitialFit(yy,Ytran+.0000,Strpidx,z0,WLofPixel,LambdaRfit,deltad,gammafit,Amplfit,bfit,Yscreenfit)
    AmplTwoXU1, AmplTwoYU1 = AmplTwoU1[0], AmplTwoU1[1]
    
    AmplTwoU2 = superamplitudeInitialFit(yy,Ytran+.0000,Strpidx,z0-deltad,WLofPixel,LambdaRfit,deltad,gammafit,Amplfit,bfit,Yscreenfit)
    phaseR = AmplTwoU2[2]
    bfitTwo = AmplTwoU2[3]
    AmplTwoXU2, AmplTwoYU2 = AmplTwoU2[0]*np.exp(1.j*2.*np.pi/WLofPixel*phaseR) + bfitTwo, AmplTwoU2[1]*np.exp(1.j*2.*np.pi/WLofPixel*phaseR) + bfitTwo

    ElectronSource = np.array(IFunc.ElectronSourcefunction(yy))

    #ElectronSource = ElectronSource/sum(ElectronSource)

    fresnelforward1 = np.array(IFunc.squaredphase(yy,Ytran,0,zpfit,WLofPixel)) #Here zero is true, it is seen as complementary to deltad + Lu
    fresnelforward2 = np.array(IFunc.squaredphase(yy,Ytran,deltad + AUndu.Lu,zpfit,WLofPixel))

    origprimeXU1 = IFunc.convolvefresnel(AmplTwoXU1,fresnelforward1,ElectronSource)/(WLofPixel*100.)
    origprimeYU1 = IFunc.convolvefresnel(AmplTwoYU1,fresnelforward1,ElectronSource)/(WLofPixel*100.)

    origprimeXU2 = IFunc.convolvefresnel(AmplTwoXU2,fresnelforward2,ElectronSource)/(WLofPixel*100.)
    origprimeYU2 = IFunc.convolvefresnel(AmplTwoYU2,fresnelforward2,ElectronSource)/(WLofPixel*100.)    

    Sourcelst1D1000 = (np.conjugate(origprimeXU1 + origprimeXU2)*(origprimeXU1 + origprimeXU2)+np.conjugate((origprimeYU1 + origprimeYU2))*(origprimeYU1 + origprimeYU2)) - np.array([.5]*240+[1.]*606+[.5]*(1000-606-240))#20220404 - np.array([.5]*240+[1.]*606+[.5]*(1000-606-240))

    Sourcelst1D1000 = Sourcelst1D1000.real
    
    Sourcelst1D1000 = Sourcelst1D1000*(1+elevationfit*deltad/100)

    #-0.005902848821248483
    #plotaperture(zetaprime2,Sourcelst1D1000)    
    return Sourcelst1D1000