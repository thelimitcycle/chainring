import numpy as np
import scipy as sp 
import finite_difference_matrix as fdm 
import matplotlib.pyplot as plt

n = 100
Sdiy = fdm.simpsons_rule_matrix(n+1)

domain = 2*np.pi
x = domain* np.arange(n+1)/n
h = domain/n

#y = np.cos(x)
y = x**2
Ydiy = 1.0/3*h*Sdiy.dot(y)

Y = sp.integrate.cumulative_simpson(y,x=x,initial=0)
Y_an = 1.0/3 * x**3

Ydiff = Ydiy-Y_an
print(Sdiy)
print(Ydiy)
print(Y)
plt.plot(x,Y,x,Ydiy)
plt.legend(['scipy', 'matrix mult'])
plt.show()

plt.figure()
plt.plot(x,Ydiff)
plt.scatter(x,Ydiff,color='red')
plt.show()

