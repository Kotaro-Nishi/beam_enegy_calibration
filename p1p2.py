from scipy.optimize import fsolve
import numpy as np

b=1.
a=2.
p1 = 0.00001
p2 = -3.

def f(phi):
    return ((b**2. - a**2.)*np.sin(phi)*np.cos(phi)+p1*a*np.sin(phi)-p2*b*np.cos(phi))
z = fsolve(f,[0.,.5,1,2,3,4,5,6,7]) 
print(z%np.pi/np.pi*180.)
print(z/np.pi*180.)

print((p1-a*np.cos(z))/(b*np.cos(z))*(b**2.*np.cos(z)**2.+a**2.*np.sin(z)**2.)**.5)
