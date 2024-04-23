import numpy as np
from scipy.sparse import diags

#https://web.media.mit.edu/~crtaylor/calculator.html

def create_stiffness_matrix(num_nodes):
    diagonals = [-1.0/12* np.ones(num_nodes-2),
                4.0/3 * np.ones(num_nodes-1),
                -5.0/2 * np.ones(num_nodes),
                4.0/3 * np.ones(num_nodes-1),
                -1.0/12* np.ones(num_nodes-2)]
    
    offsets = [-2,-1, 0, 1, 2]
    stiffness_matrix = diags(diagonals, offsets, shape=(num_nodes, num_nodes), format='csr').toarray()
    
    #boundary conditions
    stiffness_matrix[0,0:5] = 1.0/12 * np.array([35,-104,114,-56,11])
    stiffness_matrix[1,0:5] = 1.0/12 * np.array([11,-20,6,4,-1])

    #bounday conditions
    #there is a bug here!
    stiffness_matrix[num_nodes-2,-5:] = 1.0/12 * np.array([-1,4,6,-20,11])
    stiffness_matrix[num_nodes-1,-5:] = 1.0/12 * np.array([11,-56,114,-104,35])

    #I think I need some periodic boundary conditions on these matrices
    return stiffness_matrix


def create_damping_matrix(num_nodes):
    # Assuming 1D truss elements with constant material properties
    diagonals = [1.0/12* np.ones(num_nodes-2),
                -8.0/12 * np.ones(num_nodes-1),
                8.0/12 * np.ones(num_nodes-1),
                -1.0/12* np.ones(num_nodes-2)]
    offsets = [-2,-1, 1, 2]
    damping_matrix = diags(diagonals, offsets, shape=(num_nodes, num_nodes), format='csr').toarray()
    #need to think about the boundary conditions. to minimize error
    #fourth order, one side only boundary conditions
    damping_matrix[0,0]= -25.0/12
    damping_matrix[0,1]= 4
    damping_matrix[0,2] = -3
    damping_matrix[0,3] = 4.0/3
    damping_matrix[0,4] = -0.25
    #fourth order, one left, three right.... need to do the linalg on this
    damping_matrix[1,0] = -3.0/12
    damping_matrix[1,1] = -10.0/12
    damping_matrix[1,2] = 1.5
    damping_matrix[1,3] = -0.5
    damping_matrix[1,4] = 1.0/12
    #now the boundary conditions for the other side
    damping_matrix[num_nodes-1,num_nodes-5] = 0.25
    damping_matrix[num_nodes-1,num_nodes-4] = -4.0/3
    damping_matrix[num_nodes-1,num_nodes-3] = 3
    damping_matrix[num_nodes-1,num_nodes-2] = -4
    damping_matrix[num_nodes-1,num_nodes-1] = 25/12

    damping_matrix[num_nodes-2,num_nodes-5] = -1.0/12
    damping_matrix[num_nodes-2,num_nodes-4] = 0.5
    damping_matrix[num_nodes-2,num_nodes-3] = -1.5
    damping_matrix[num_nodes-2,num_nodes-2] = 5.0/6
    damping_matrix[num_nodes-2,num_nodes-1] = 1.0/4

    #print(damping_matrix.shape)
    return damping_matrix