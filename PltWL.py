import matplotlib.pyplot as plt
import numpy as np

from scipy.optimize import curve_fit

z0 = 10.
z1 = 3.03
z2 = .32
zp = z1 + z2
zc = .68
zcp = .4

za = z0 + zp
zb = zc + zcp

print(za**2 * zcp / z0 * (1/(z0*zp*(zc+zcp-zp)-zp**2*zc-zb**2*z0)))
print(1/(za**2 * zcp / z0 * (1/(z0*zp*(zc+zcp-zp)-0.*zp**2*zc-zb**2*z0))))
#exit()

WLths = [3.983459533328215e-07, 3.9851425990667447e-07, 3.9868256648052747e-07, 3.988508730543805e-07, 3.990191796282335e-07, 3.991874862020865e-07, 3.993557927759395e-07, 3.995240993497925e-07, 3.996924059236455e-07, 3.9986071249749847e-07, 4.0002901907135147e-07, 4.0019732564520446e-07, 4.0036563221905745e-07, 4.0053393879291045e-07, 4.007022453667635e-07, 4.008705519406165e-07, 4.010388585144695e-07, 4.0120716508832247e-07, 4.0137547166217546e-07, 4.0154377823602846e-07, 4.0171208480988145e-07, 4.0188039138373444e-07, 4.0204869795758744e-07, 4.0221700453144043e-07, 4.023853111052934e-07, 4.0255361767914647e-07, 4.0272192425299946e-07, 4.0289023082685246e-07, 4.0305853740070545e-07, 4.0322684397455844e-07, 4.0339515054841144e-07, 4.0356345712226443e-07, 4.037317636961174e-07, 4.039000702699704e-07, 4.040683768438234e-07, 4.042366834176764e-07, 4.0440498999152945e-07, 4.0457329656538244e-07, 4.0474160313923544e-07, 4.0490990971308843e-07, 4.050782162869414e-07, 4.052465228607944e-07, 4.054148294346474e-07, 4.055831360085004e-07, 4.057514425823534e-07, 4.059197491562064e-07, 4.060880557300594e-07, 4.0625636230391243e-07, 4.064246688777654e-07, 4.065929754516184e-07, 4.067612820254714e-07, 4.069295885993244e-07, 4.070978951731774e-07, 4.072662017470304e-07, 4.074345083208834e-07, 4.076028148947364e-07, 4.0777112146858937e-07, 4.0793942804244236e-07, 4.081077346162954e-07, 4.082760411901484e-07, 4.084443477640014e-07, 4.086126543378544e-07, 4.087809609117074e-07, 4.0894926748556037e-07, 4.0911757405941337e-07]

gammas = [
352.4449542547033,
352.42027829085185,
352.432152150499,
352.43895144013345,
352.4289591098259,
352.43846279453146,
352.4399427973816,
352.4394134427969,
352.44653807947645,
352.4371638089255,
352.4483363525002,
352.44653075982154,
352.4463699371771,
352.45599326173874,
352.4503209141664,
352.4542910525692,
352.45388689042386,
352.4517798427986,
352.4504699098189,
352.43823473306054,
352.4545308244095,
352.4541134781777,
352.4500783260651,
352.44535823283707,
352.447851945601,
352.42380070947917,
352.45045743814364,
352.4324488111239,
352.45299887629625,
352.43235766792793,
352.4560833380098,
352.43885175702036,
352.4280659569531,
352.4252903461915,
352.44193207099784,
352.42554385285575,
352.43520733906615,
352.41361234767083,
352.42431748192973,
352.43318936445974,
352.45207496706155,
352.43000197427546,
352.4333427330366,
352.4681871923539,
352.46119659802895,
352.4474679119111,
352.4532388375156,
352.45114625152485,
352.45363455357517,
352.43908889829044,
352.46042786374943,
352.48767967814405,
352.49355886040485,
352.488834817007,
352.51467755497384,
352.53484449568606,
352.51832064723334,
352.5094570603793,
352.5162039370611,
352.5525803849457,
352.54396429604293,
352.5547547248807,
352.5626176156441,
352.5514328257577,
352.56316892082197
]


def gaussian(x,sigma,mu):
    return 1./(sigma*np.sqrt(2*np.pi)) *np.exp(-(x-mu)**2./(2*sigma**2.))

def linfunc(x,m,b):
    return m*x+b

gammas_good = gammas[0:len(gammas)-14]
WLths_good = WLths[0:len(WLths)-14]

gammas_bad = gammas[len(gammas)-14:-1]
WLths_bad = WLths[len(WLths)-14:-1]

num_bins = 18
fntsze = 30
fig, ax = plt.subplots(figsize =(10, 10))
# Creating plot
n, bins, patches = plt.hist(gammas_good, num_bins, facecolor = 'green')
plt.title("Histogram of gamma over the spectrum", fontsize=35)

plt.xticks(fontsize=30)
plt.yticks(fontsize=30)

m = np.mean(gammas_good)
sig = np.std(gammas_good)

binlin = np.linspace(bins[-1],bins[0],100)

#ax.set_title('Sales by City', fontsize=15, color= 'blue', fontweight='bold')
ax.set_xlabel('gamma',fontsize = fntsze)
ax.set_ylabel('counts',fontsize = fntsze)
ax.plot(binlin,gaussian(np.array(binlin),sig,m), linewidth = 4)

  
# show plot
plt.tight_layout() 
plt.show()
#exit()
#print(np.std(np.array(gammas_good)))
#0.004808759252942803
print((WLths_good[-1]-WLths_good[0])/0.004808759252942803)

WLths_good = np.array(WLths_good)
WLths_bad = np.array(WLths_bad)

lowWL = 404.6565
highWL = 407.8737

#lowg = 352.4
#highg = 352.6

lowg = 352.3
highg = 352.5


print(WLths_good*10.e8)


#fntsze = 30
fig = plt.figure()
#plt.title("gamma vs. wavelength", fontsize=35)

ax = fig.add_subplot(111)

label_size = 30
ax.tick_params(labelsize=label_size)
#ax.set_xlim([402.5,408])
#ax.set_xlim([402.5,408])
#ax.set_ylim([lowg,highg])
ax.plot((np.array(WLths)*10.e8)[0:51],gammas[0:51],'b.', markersize=15)
ax.plot([lowWL,lowWL],[lowg,highg],'r')
#ax.plot([highWL,highWL],[lowg,highg],'r')
ax.set_xlabel('Wavelength [nm]',fontsize = fntsze)
ax.set_ylabel('gamma',fontsize = fntsze)
plt.show()


#fntsze = 30
fig = plt.figure()
#plt.title("gamma vs. wavelength", fontsize=35)

ax = fig.add_subplot(111)

label_size = 30
ax.tick_params(labelsize=label_size)
#ax.set_xlim([402.5,408])
#ax.set_xlim([402.5,408])
#ax.set_ylim([lowg,highg])
ax.plot(WLths_good*10.e8,gammas_good,'b.', markersize=15)
ax.plot(WLths_bad*10.e8,gammas_bad,'r.', markersize=15)
ax.plot([lowWL,lowWL],[lowg,highg],'r')
ax.plot([highWL,highWL],[lowg,highg],'r')
ax.set_xlabel('Wavelength [nm]',fontsize = fntsze)
ax.set_ylabel('gamma',fontsize = fntsze)
plt.show()


guess = np.array([0.,350.])


lowg = 352.
highg = 353.



popt,pcov = curve_fit(linfunc,np.array(WLths_good)*10.e8-WLths_good[0]*10.e8,np.array(gammas_good)*.511,guess*.511)
#xargs = np.linspace()
print(popt)
print(np.sqrt(pcov[0][0]))
print(np.sqrt(pcov[1][1]))

fig = plt.figure()
#plt.title("gamma vs. wavelength", fontsize=35)
ax = fig.add_subplot(111)

label_size = 30
ax.tick_params(labelsize=label_size)
#ax.set_xlim([402.5,408])
#ax.set_xlim([402.5,408])
ax.set_ylim([lowg*.511,highg*.511])
ax.plot(WLths_good*10.e8,np.array(gammas_good)*.511,'b.', markersize=15)
#ax.plot(WLths_bad*10.e8,gammas_bad,'r', markersize=15)
ax.plot(np.array(WLths_good)*10.e8,linfunc(np.array(WLths_good*10.e8-WLths_good[0]*10.e8),*popt),'r', linewidth=4)

ax.plot([lowWL,lowWL],[(lowg+.2)*.511,(highg-.2)*.511],'k-.',linewidth=5)
#ax.plot([highWL,highWL],[lowg,highg],'r')
ax.set_xlabel('Wavelength [nm]',fontsize = fntsze)
ax.set_ylabel('Energy [MeV]',fontsize = fntsze)
ax.legend(['Data','Linear Fit','Calibration Wavelength'])
plt.show()

exit()
popt,pcov = curve_fit(linfunc,np.array(WLths_good)*10.e8-WLths_good[0]*10.e8,np.array(gammas_good),guess)
#xargs = np.linspace()
print(popt)
print(np.sqrt(pcov[0][0]))
print(np.sqrt(pcov[1][1]))

fig = plt.figure()
#plt.title("gamma vs. wavelength", fontsize=35)
ax = fig.add_subplot(111)

label_size = 30
ax.tick_params(labelsize=label_size)
#ax.set_xlim([402.5,408])
#ax.set_xlim([402.5,408])
ax.set_ylim([lowg,highg])
ax.plot(WLths_good*10.e8,gammas_good,'b.', markersize=15)
#ax.plot(WLths_bad*10.e8,gammas_bad,'r', markersize=15)
ax.plot(np.array(WLths_good)*10.e8,linfunc(np.array(WLths_good*10.e8-WLths_good[0]*10.e8),*popt),'r', linewidth=4)

ax.plot([lowWL,lowWL],[lowg+.2,highg-.2],'k-.',linewidth=5)
#ax.plot([highWL,highWL],[lowg,highg],'r')
ax.set_xlabel('Wavelength [nm]',fontsize = fntsze)
ax.set_ylabel('gamma',fontsize = fntsze)
ax.legend(['Data','Linear Fit','Calibration Wavelength'])
plt.show()
