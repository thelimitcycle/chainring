import numpy as np
from finite_difference_matrix import stiffness_matrix
from finite_difference_matrix import damping_matrix
from finite_difference_matrix import simpsons_rule_matrix
import matplotlib.pyplot as plt

n = 100
domain = 2*np.pi
x = domain* np.arange(n)/n
y = np.cos(x)


C = damping_matrix(n)
print(C)
K = stiffness_matrix(n)
S = simpsons_rule_matrix(n)
h = domain * 1/n

yp = 1.0/h * C.dot(y)
ypp = 1.0/(h*h)*K.dot(y)
Y = 1.0/3 *h *S.dot(y)
Y_an = np.sin(x)
Yerr = Y-Y_an

plt.figure()
plt.plot(x,y,x,yp,x,ypp,x,Y)
plt.ylabel('y')
plt.xlabel('x')
plt.legend(['y', 'dy/dx', 'd2y/dx2', 'Y'])
plt.title('function and its first and second derivatives and antiderivative')
plt.show()

plt.figure()
plt.plot(x,Yerr)
plt.title("error between the numeric and analytic integral")
plt.show()