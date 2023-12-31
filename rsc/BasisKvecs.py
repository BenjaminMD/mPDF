import numpy as np

def get_Tb1():
    mu = 0.29
    k13 = np.array([1/2, mu, 0])

    Tb1_Atom1 = np.array([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ])

    Tb1_Atom2 = np.array([
        [1, 0, 0],
        [0, -1, 0],
        [0, 0, 1],
        
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        
        [-1, 0, 0],
        [0, 1, 0],
        [0, 0, -1],
    ]) * np.exp(1j*np.pi*mu)


    Tb1_Atom3 = np.array([
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        
        [1, 0, 0],
        [0, -1, 0],
        [0, 0, -1],
        
        [-1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],    
        
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],

    ]) * np.exp(1j*np.pi*(mu-0.5))

    Tb1_Atom4 = np.array([
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        
        [-1, 0, 0],
        [0, -1, 0],
        [0, 0, 1],
        
        [-1, 0, 0],
        [0, -1, 0],
        [0, 0, 1],    
        
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],

    ]) * np.exp(1j*np.pi/2)
    
    return k13, Tb1_Atom1, Tb1_Atom2, Tb1_Atom3, Tb1_Atom4
