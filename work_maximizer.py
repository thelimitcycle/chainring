import numpy as np
from scipy.sparse import diags
from scipy.integrate import simpson
import matplotlib.pyplot as plt
import pandas as pd
import os
from finite_difference_matrix import damping_matrix 
from finite_difference_matrix import stiffness_matrix
from finite_difference_matrix import simpsons_rule_matrix

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
    fn = force_data['Fn'].to_numpy()[:-1]
    fn = fn + np.roll(fn,180)
    r = 0.10
    r_regular = np.ones(len(fn))*r

    num_angles = len(fn)
    C_ang = damping_matrix(num_angles)
    K_ang = stiffness_matrix(num_angles)

    #print(fn)
    theta = 2*np.pi*np.arange(360)/360
    h = 2*np.pi / num_angles

    #print(y.shape)
    plt.figure()

    fnp = (1.0/h)*C_ang.dot(fn)
    arc_length_force = np.sqrt(fn*fn + fnp*fnp)
    #print(arc_length_force)
    round_perimeter = 2*np.pi*r
    l = 1/(2*round_perimeter) * simpson(arc_length_force,x=theta)
    
    denom = 2*l
    average_force = np.mean(fn)
    print("average force: " + str(average_force))
    print("lagrange multiplier: " + str(denom))
    #print(fn)
    #print(average_force)
    r_optimized = r*fn/average_force
    torque_original = fn*r_regular
    torque_optimized = fn*r_optimized

    work_original = np.abs(simpson(torque_original,x=theta))
    work_optimized = np.abs(simpson(torque_optimized,x=theta))

    print("Original work = " + str(work_original))
    print("optimized work = " + str(work_optimized))
    rpo = C_ang.dot(r_optimized)
    #print(str(r_optimized.max()))
    #print(str(r_optimized.min()))
    print("circumference of original gear (theoretical): " + str(round_perimeter))
    print("circumference of original gear (numerical): " + str(simpson(r_regular,x=theta)))
    print("circumference of optimized gear (numerical): " + str(simpson(np.sqrt(r_optimized*r_optimized + rpo*rpo),x=theta)))
    #print(yp)
    #print(yp.shape)
    #plt.plot(theta,r_regular, theta,r_optimized)
    #plt.figure()

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
    plt.scatter(x_opt,y_opt)
    plt.title('chainring shape for '+ profile_name)
    plt.ylim([-.15,.15])
    plt.xlim([-.15,.15])
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    plt.savefig('plots/perimeter/' + profile_name + "_perimeter.png")
    #plt.show()

    results_dict ={'profile':profile_name,
                   'round_work':work_original,
                   'optimized_work': work_optimized}
    return results_dict

if __name__ ==  "__main__":
    fnames = os.listdir("kautz/")
    col_names = ['profile', 'round_work','optimized_work']
    df = pd.DataFrame(columns=col_names)
    for fname in fnames:
        result = work_optimizer("kautz/" + fname)
        df.loc[len(df.index)] = result

    df['percentage_difference'] = 100*(df['optimized_work'] - df['round_work']).div(df['round_work'])
    print(df)