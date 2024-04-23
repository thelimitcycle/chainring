import numpy as np
from finite_difference_matrix import create_stiffness_matrix
from finite_difference_matrix import create_damping_matrix
import matplotlib.pyplot as plt

#so these are experiments with the matrices and with calculus in general in polar coordinates
# I think I'm going to have to add some differential constraints on the Euler-LaGrange equation
# so I'm practicing the numerics here

#get my angle
n = 360
theta = 2* np.pi*np.arange(n)/n 

#just a circle
r = 1*np.ones(len(theta))


fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(theta, r)
ax.set_rmax(2)
ax.set_rticks([0.5, 1, 1.5, 2])  # Less radial ticks
ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
ax.grid(True)

ax.set_title("A line plot on a polar axis", va='bottom')
plt.show()

#lets plot something more exagerated
#https://math.stackexchange.com/questions/315386/ellipse-in-polar-coordinates
a = 2
b = 1

r2 = a*b / np.sqrt((b*np.cos(theta))**2 + (a*np.sin(theta))**2) # oval
#r2 = np.sin(2*theta)

#lets see if the dr/dtheta looks like anything...

#I think I need som periodic bounday conditions on the matrices
C = create_damping_matrix(n)
h = 2*np.pi / n
r2p = (1.0/h)*C.dot(r2)
r2pp2 = (1.0/h)*C.dot(r2p)
K = create_stiffness_matrix(n)

#wait f, I think I need to do the chain rule
#https://mechref.engr.illinois.edu/dyn/rvy.html maybe

#I think there is a bug here or boundary conditions aren't correct
#https://www.cs.auckland.ac.nz/courses/compsci369s1c/lectures/GG-notes/Z-archive2010/CS369-FDM-FEM.pdf
# I think I need to implement 2nd order 
r2pp = (1/h**2)*K.dot(r2)

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(theta, r2,theta, r2p,theta, r2pp)
ax.set_rmax(2)
ax.set_rticks([0.5, 1, 1.5, 2])  # Less radial ticks
ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
ax.grid(True)

ax.set_title("A line plot on a polar axis", va='bottom')

plt.show()
plt.figure()
plt.plot(theta,r2, theta,r2p,theta,r2pp)
plt.show()

#ok now this should be fun. 