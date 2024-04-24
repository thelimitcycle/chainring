import numpy as np
from finite_difference_matrix import create_stiffness_matrix
from finite_difference_matrix import create_damping_matrix
import matplotlib.pyplot as plt

n = 100
domain = 4
x = domain* np.arange(n)/n
y = x**2

C = create_damping_matrix(n)
print(C)
K = create_stiffness_matrix(n)
h = domain * 1/n

yp = 1.0/h * C.dot(y)
ypp = 1.0/(h*h)*K.dot(y)

plt.figure()
plt.plot(x,y,x,yp,x,ypp)
plt.ylabel('y')
plt.xlabel('x')
plt.legend(['y', 'dy/dx', 'd2y/dx2'])
plt.title('function and its first and second derivatives')
plt.show()