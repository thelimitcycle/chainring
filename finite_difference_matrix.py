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


def simpsons_rule_matrix(num_nodes):
    #I think its a lower triangular matrix with the columns being the simpson's coefficients
    # wait I'm not sure if the intragral as a function will work. 
    #simpson's rule matrix goes here
    # a matrix approach to numerical intrgration maybe
    # thats right the fancy name is 'quadrature'

    # may have to implement a numerical quadrature algorithm and use polynomial or cosine 
    # basis functions. like semhat
    #https://www3.nd.edu/~coast/jjwteach/www/www/60130/New%20Lecture%20Notes_PDF/CE60130_Lecture%2013%20with_footer-v04.pdf


    #I think I can do it with simpson's rule
    n = num_nodes

    S = np.ones((n,n))
    for k in range(n):
        if k == 0:
            pass
        elif k == n-1:
            pass
        elif k%2 == 1:
            S[:,k] = S[:,k]*4
        else:
            S[:,k] = 2*S[:,k]
        S[k,k] = 1
    S = np.tril(S)

    for k in range(n):
        if k%2 != 0:
            S[k,:] = np.roll(S[k-1,:],1)

    #print(S)
    return(S)