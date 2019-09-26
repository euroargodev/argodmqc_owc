import numpy as np
from scipy import interpolate

mymin,mymax = 0,3
X = np.linspace(mymin,mymax,4)
Y = np.linspace(mymin,mymax,4)

x,y = np.meshgrid(X,Y)

test = np.array([[ 1.2514318 ,  1.25145821,  1.25148472,  1.25151133],
       [ 1.25087456,  1.25090105,  1.25092764,  1.25095435],
       [ 1.25031581,  1.25034238,  1.25036907,  1.25039586],
       [ 1.24975557,  1.24978222,  1.24980898,  1.24983587]])

f = interpolate.interp2d(x,y,test,kind='cubic')

# use linspace so your new range also goes from 0 to 3, with 8 intervals
Xnew = np.linspace(mymin,mymax,8)
Ynew = np.linspace(mymin,mymax,8)

test8x8 = f(Xnew,Ynew)

print(test8x8)