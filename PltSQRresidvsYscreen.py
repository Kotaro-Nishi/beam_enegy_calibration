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

gammas = [
484027.88,
483554.84,
467514.1,
445445.28,
452850.6,
442865.03,
439107.34,
441357.97,
409249.94,
404856.16,
406952.22,
466836.25,
485799.06,
361772.4,
351793.28,
336979.66,
338373.5,
335179.62,
325658.44,
318380.8,
328826.34,
333406.22,
331089.88,
338856.8,
334477.53,
319284.84,
304767.7,
300072.7,
308095.8,
315395.7,
314854.44,
335731.7,
348978.66,
347914.75,
353210.1,
414333.38,
441050.6,
369867.2,
363201.12,
365506.38,
365389.22,
377821.47,
382310.72,
407299.53,
418110.5
]

gammastds = [
3960.368,
4395.057,
5029.6978,
3652.0076,
7021.3677,
6101.7075,
6685.5996,
4533.62,
21156.367,
3312.214,
27096.494,
4468.851,
33852.125,
2287.9502,
1731.0444,
3337.716,
2092.9138,
3190.8687,
3019.4553,
1994.2175,
2290.4592,
2311.2185,
2024.9574,
2365.472,
2292.6692,
2019.6035,
1740.0579,
1796.405,
2783.989,
3221.6572,
2100.1235,
3701.8328,
6304.723,
2621.4275,
4964.822,
17693.828,
24132.527,
28667.52,
1590.9558,
2304.9077,
10637.269,
15547.184,
8248.665,
3480.8784,
3141.0002
]

WLths = np.array(range(len(gammas)))*10+190.

fntsze = 20

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title("Sum of squared residuals vs. Position on CMOS screen", fontsize=25)
#ax.title("Sum of squared residuals vs. Position on CMOS screen")#, fontsize=35)
ax.tick_params(axis='x', labelsize=15)
ax.tick_params(axis='y', labelsize=15)
#ax.xticks(fontsize=30)
#ax.yticks(fontsize=30)

#ax.set_xlim([402.5,408])
#ax.set_xlim([402.5,408])
#ax.set_ylim([25E4,55E4])
ax.set_xlabel(r'Position on CMOS screen [$\mu$m]',fontsize = fntsze)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax.set_ylabel('Sum of squared residuals',fontsize = fntsze)
#ax.fill_between(WLths, gammas-10.*np.array(gammastds), gammas+10.*np.array(gammastds))
ax.plot(WLths,gammas,'r')#, markersize=15)

plt.show()

exit()

gammas_good = gammas[0:len(gammas)-40]
WLths_good = WLths[0:len(WLths)-40]

gammas_bad = gammas[len(gammas)-40:-1]
WLths_bad = WLths[len(WLths)-40:-1]

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
plt.title("gamma vs. wavelength", fontsize=35)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
ax = fig.add_subplot(111)
#ax.set_xlim([402.5,408])
ax.set_xlim([402.5,408])
ax.set_ylim([lowg,highg])
ax.plot(WLths_good*10.e8,gammas_good,'b.', markersize=15)
ax.plot(WLths_bad*10.e8,gammas_bad,'r.', markersize=15)
ax.plot([lowWL,lowWL],[lowg,highg],'r')
ax.plot([highWL,highWL],[lowg,highg],'r')
ax.set_xlabel('Wavelength [nm]',fontsize = fntsze)
ax.set_ylabel('gamma',fontsize = fntsze)
plt.show()


guess = np.array([0.,350.])


lowg = 350.
highg = 354.


popt,pcov = curve_fit(linfunc,np.array(WLths_good)*10.e8-WLths_good[0]*10.e8,np.array(gammas_good),guess)
#xargs = np.linspace()
print(popt)
print(np.sqrt(pcov[0][0]))
print(np.sqrt(pcov[1][1]))

fig = plt.figure()
plt.title("gamma vs. wavelength", fontsize=35)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
ax = fig.add_subplot(111)
#ax.set_xlim([402.5,408])
ax.set_xlim([402.5,408])
ax.set_ylim([lowg,highg])
ax.plot(WLths_good*10.e8,gammas_good,'b', markersize=15)
ax.plot(WLths_bad*10.e8,gammas_bad,'r', markersize=15)
ax.plot(np.array(WLths_good)*10.e8,linfunc(np.array(WLths_good*10.e8-WLths_good[0]*10.e8),*popt),'r', markersize=15)

ax.plot([lowWL,lowWL],[lowg,highg],'r')
ax.plot([highWL,highWL],[lowg,highg],'r')
ax.set_xlabel('Wavelength [nm]',fontsize = fntsze)
ax.set_ylabel('gamma',fontsize = fntsze)
plt.show()
