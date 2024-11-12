import numpy as np
from scipy.sparse import diags
from scipy.integrate import simpson
import matplotlib.pyplot as plt
import pandas as pd
import os
from finite_difference_matrix import damping_matrix 
from finite_difference_matrix import stiffness_matrix
from finite_difference_matrix import simpsons_rule_matrix


# function for taking the curvature of a function. 
# r is a function of theta and I think it should be equally sampled
# yeah it should because of the finite differences
# make sure it lives [0, 2*pi]
def get_curvature(r, C, K):
    n = len(r)
    #make our first and second order differentiation matrices
    # nah take these as parameters
    # C = damping_matrix(n)
    # K = stiffness_matrix(n)
    
    #grid spacing
    h = 2*np.pi / n
    #take derivates of r(theta) with respect to theta
    drdtheta = (1.0/h)*C.dot(r)
    d2rdtheta2 = (1.0/h**2)*K.dot(r)
    #formula for curvature. first with real numbers
    #https://mathworld.wolfram.com/Curvature.html
    #http://mathonline.wikidot.com/the-curvature-of-plane-polar-curves
    # ok so this formula doesn't really work well on sampled data.  :/ 
    # maybe need to smooth the data out somehow.
    curvy = (2*drdtheta*drdtheta - r*d2rdtheta2 + r*r) / np.power((drdtheta*drdtheta + r*r),1.5)

    #I'm, not sure how to smoothn out the curvature signal
    #
    #smooth out the curvature signal.
    #https://dsp.stackexchange.com/questions/63885/real-time-numerical-differentiation-of-signals
    #
    return curvy



# Example parameters
def work_optimizer(f_path):
    #TODO:
    #1. create a dataframe for keeping track of the following parameters
    #   -length
    #   -work
    #   -percentage difference from round
    #2. make plots of each of the gear profiles and the efficiency changes.
    #3. i think the best way to do this is to return a dictionary and then append the dictionary to a dataframe
    #in the main function
    profile_name = f_path.split(r'/')[1][:-4]

    force_data = pd.read_csv(f_path,skiprows=[0,1,2,3,4,5,6,8],delim_whitespace=True)
    # normal component of pedal force
    fn = force_data['Fn'].to_numpy()[:-1]
    # the dataset measures one foot pedaling
    # assuming symetry, offset the force curve by pi radians for both feet
    fn = fn + np.roll(fn,180)
    
    #radius in meters
    r = 0.10

    #circular gear
    r_regular = np.ones(len(fn))*r

    num_angles = len(fn)
    C_ang = damping_matrix(num_angles)
    K_ang = stiffness_matrix(num_angles)

    #print(fn)
    theta = 2*np.pi*np.arange(360)/360
    h = 2*np.pi / num_angles

    fnp = (1.0/h)*C_ang.dot(fn)
    #this is the integrand of the arc length formula in polar coordinates
    arc_length_force = np.sqrt(fn*fn + fnp*fnp)
    #print(arc_length_force)
    round_perimeter = 2*np.pi*r

    l = 1/(2*round_perimeter) * simpson(arc_length_force,x=theta)
    
    denom = 2*l # i guess this is the consant
    average_force = np.mean(fn)
    print("average force: " + str(average_force))
    print("lagrange multiplier: " + str(denom))
    #print(fn)

    #the solution to the differential equation is acutally really simple lol
    # normalize for average force. intermediate value theorem or something
    # i think this is a dimensional units argument. Idk if it's valid though
    # need to step through the algebra
    r_optimized = r*fn/average_force 
    torque_original = fn*r_regular
    torque_optimized = fn*r_optimized

    work_original = np.abs(simpson(torque_original,x=theta))
    work_optimized = np.abs(simpson(torque_optimized,x=theta))

    #now need to check for negative curvature
    curvy = get_curvature(r_optimized, C_ang, K_ang)
    print("length of curvature is "  + str(len(curvy)) + "length of r_optimized is " + str(len(r_optimized)))
    negative_curvature = curvy > 0


    print("Original work = " + str(work_original))
    print("optimized work = " + str(work_optimized))
    rpo = C_ang.dot(r_optimized)
    #print(str(r_optimized.max()))
    #print(str(r_optimized.min()))
    print("circumference of original gear (theoretical): " + str(round_perimeter))
    print("circumference of original gear (numerical): " + str(simpson(r_regular,x=theta)))
    print("circumference of optimized gear (numerical): " + str(simpson(np.sqrt(r_optimized*r_optimized + rpo*rpo),x=theta)))

    plt.figure()
    plt.plot(theta,torque_original,theta,torque_optimized)
    plt.legend(["torque with round gear", "torque with optimized gear"])
    plt.xlabel("crank angle (rad)")
    plt.ylabel("torque (Nm)")
    plt.xlim([0, 2*np.pi])
    plt.ylim([-75,0])
    plt.title('Original and optimized torque values for ' + profile_name)
    plt.savefig('plots/torque/' + profile_name +'_torque.png')
    #plt.show()

    #eeeh I'll just transform to cartesian coordinates and then do a scatter
    x_opt = r_optimized*np.sin(theta)
    y_opt = r_optimized*np.cos(theta)
    
    plt.figure(figsize=(8,8))
    plt.scatter(x_opt,y_opt,c=negative_curvature)
    plt.title('chainring shape for '+ profile_name)
    plt.ylim([-.15,.15])
    plt.xlim([-.15,.15])
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    plt.savefig('plots/perimeter/' + profile_name + "_perimeter.png")
    #plt.show()

    #plot the curvature
    # Create a plot with two y-axes
    fig, ax1 = plt.subplots()

    # Plot on the first y-axis
    ax1.plot(theta, r_optimized, 'b-', label='sin(x)')
    ax1.set_xlabel('r_optimized')
    # ax1.set_ylabel('sin(x)', color='b')
    ax1.tick_params(axis='y', labelcolor='b')

    # Create a second y-axis sharing the same x-axis
    ax2 = ax1.twinx()
    ax2.plot(theta, curvy, 'r-', label='cos(x)')
    # ax2.set_ylabel('cos(x)', color='r')
    ax2.tick_params(axis='y', labelcolor='r')

    # Show plot
    fig.tight_layout()  # Adjust layout to avoid overlap
    plt.title("r_optimized and curvature")
    plt.savefig('plots/curvature/' + profile_name + "_curvature.png")

    results_dict ={'profile':profile_name,
                   'round_work':work_original,
                   'optimized_work': work_optimized,}
    
    return [results_dict, r_optimized]

if __name__ ==  "__main__":
    fnames = os.listdir("kautz/")
    col_names = ['profile', 'round_work','optimized_work']
    df = pd.DataFrame(columns=col_names)
    df2 = pd.DataFrame(columns=fnames)
    for fname in fnames:
        result = work_optimizer("kautz/" + fname)
        df.loc[len(df.index)] = result[0]
        df2[fname] = result[1]

    df['percentage_difference'] = 100*(df['optimized_work'] - df['round_work']).div(df['round_work'])
    print(df)

    df2.to_csv('r_optmized.csv')