import numpy as np
from finite_difference_matrix import stiffness_matrix
from finite_difference_matrix import damping_matrix
import matplotlib.pyplot as plt

#so these are experiments with the matrices and with calculus in general in polar coordinates
# I think I'm going to have to add some differential constraints on the Euler-LaGrange equation
# so I'm practicing the numerics here

#get my angle

#ffff I'm going to have to set up a virtual environment
#https://stackoverflow.com/questions/75602063/pip-install-r-requirements-txt-is-failing-this-environment-is-externally-mana/75696359#75696359
n = 3600
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
#plot an ellipse
a = 2
b = 1

#playing with different values
r2 = np.sin(theta)
#r2 = a*b / np.sqrt((b*np.cos(theta))**2 + (a*np.sin(theta))**2) # oval
#r2 = np.sin(2*theta)

#lets see if the dr/dtheta looks like anything...

#I think I need som periodic bounday conditions on the matrices
C = damping_matrix(n)
K = stiffness_matrix(n)

#grid spacing
h = 2*np.pi / n

#take derivates of r(theta) with respect to theta
drdtheta = (1.0/h)*C.dot(r2)
d2rdtheta2 = (1.0/h**2)*K.dot(r2)

# I think I need to implement 2nd order 

#polar derivatives:
#https://www.whitman.edu/mathematics/calculus_online/section10.02.html
#its a good resource. I'm going to program this prof's formulae for first and second derivatives

# the C and K matrices should be for the derivatives with respect to theta. 

dydx = (r2*np.cos(theta) + drdtheta*np.sin(theta))/(-1.0*r2*np.sin(theta) + drdtheta*np.cos(theta))
#I think this guy's derivation is incomplete for the second derivative
#he says its (d(dydx)/d(theta)) / (dx/d(theta))
#but x = r*cos(theta) 
#so dx/d(theta) = dr/dtheta * cos(theta) - r(theta)*sin(theta)

#using numerics in the numerator
d2ydx2 = (1.0/h)*C.dot(dydx)/(drdtheta*np.cos(theta) - r2*np.sin(theta))

#can this all be one giant matrix?
## this my idea. I know I can take the derivative with a matrix
## I should be able to integrate with a matrix. I can take the second derivative with a matrix
## so I should be able to have some sort of over-contrained set of linear equations 
## where I try and use a linear least square fit to a series of difference equations for the 
## equiperimetric contraint, and the concavity constraints. 

## it should be possible I just need to setup the matrix. It's going to be a weird one.

#lets see if I can find the points where the second derivative changes sign. 
# is there anything else that I can use?

#maybe gradient descent?
#https://en.wikipedia.org/wiki/Gradient_descent


#formula for curvature. first with real numbers
#https://mathworld.wolfram.com/Curvature.html
#http://mathonline.wikidot.com/the-curvature-of-plane-polar-curves
curvature = (2*drdtheta*drdtheta - r*d2rdtheta2 + r*r) / np.power((drdtheta*drdtheta + r*r),1.5)

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

ax.plot(theta, r2)
ax.set_rmax(2)
ax.set_rticks([0.5, 1, 1.5, 2])  # fewer radial ticks
ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
ax.grid(True)

ax.set_title("an oval", va='bottom')

#plot in rectangular coordinates. 

#looks like the second derivative is showing some sin(1/x) like behavior where it is unstable
plt.show()

y = r2*np.sin(theta)
x = r2*np.cos(theta)

plt.figure()
plt.plot(theta,r2, theta,dydx, theta, d2ydx2, theta, curvature)
plt.ylim([-4,4])
plt.legend(['radius (r)', 'dy/dx instantaneous slope', 'd2y/dx2 second derivative', 'curvature'])
plt.title('an oval and its first and second derivatives')
plt.xlabel('angle [rad]')
plt.ylabel('r')
plt.show()

#ok now this should be fun. 
