import numpy as np 

def get_rotation_matrix():
    a = 2*np.random.rand(3)-1
    a /= np.linalg.norm(a)
    b = np.cross(a, 2*np.random.rand(3)-1)
    b /= np.linalg.norm(b)
    c = np.cross(a, b)
    c /= np.linalg.norm(c)

    return np.transpose(np.array([a, b, c]))