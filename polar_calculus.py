import numpy as np
from finite_difference_matrix import create_stiffness_matrix
from finite_difference_matrix import create_damping_matrix
import matplotlib.pyplot as plt

#so these are experiments with the matrices and with calculus in general in polar coordinates
# I think I'm going to have to add some differential constraints on the Euler-LaGrange equation
# so I'm practicing the numerics here

#get my angle
theta = 2* np.pi*np.arange(360)/360 

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

