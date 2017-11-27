from scipy.interpolate import RegularGridInterpolator as RGI
import numpy as np 

def project_fst_ft(mol):
	'''
	mol is an NXNXN array that contains the electron density of the protein thing
	R is the rotation matrix that symbolizes the viewing direction of the microscope
	'''
	
	#print('Now running 3D fourier Transform ...')
	N = len(mol)
	print ('okkk')
	rho_hat = np.fft.fftshift(np.fft.fftn(mol))
	wx, wy, wz = np.meshgrid(np.arange(-(N-1)/2, (N-1)/2+1, dtype= int), np.arange(-(N-1)/2, (N-1)/2+1, dtype= int), np.arange(-(N-1)/2, (N-1)/2+1, dtype= int))
	rho_hat = ((-1)**np.abs(wx+wy+wz)/(N**3)*rho_hat)

	#print('Now running Linear Interpolation ...')
	N_range = np.arange(-(N-1)/2, (N-1)/2+1)
	rho_hat = RGI((N_range,N_range,N_range), rho_hat, bounds_error= False, fill_value=0)

	### make 2D mesh that takes slice of rho_hat which NxNx3 grid 
	#print('Creating 2D Sampling grid ...')
	eta_x, eta_y = np.meshgrid(np.arange(-(N-1)/2, (N-1)/2+1, dtype= int), np.arange(-(N-1)/2, (N-1)/2+1, dtype= int), indexing= "ij")
	eta_x= eta_x[...,np.newaxis]
	eta_y= eta_y[...,np.newaxis]

	return (eta_x, eta_y, )

def project_fst_rotation(R, mol_ft)
	R= np.transpose(R)
	a= R[0]
	b= R[1]
	grid = eta_x*a + eta_y*b
	sample_grid = rho_hat(grid)

	em_slice = (-1)**np.abs(eta_x[...,0]+eta_y[...,0])*(N**2)*sample_grid

	em_slice = np.fft.ifftn(np.fft.ifftshift(em_slice))

	return np.real(em_slice)



def project_fst(R, mol):
	'''
	mol is an NXNXN array that contains the electron density of the protein thing
	R is the rotation matrix that symbolizes the viewing direction of the microscope
	'''
	
	#print('Now running 3D fourier Transform ...')
	N = len(mol)
	print ('okkk')
	rho_hat = np.fft.fftshift(np.fft.fftn(mol))
	wx, wy, wz = np.meshgrid(np.arange(-(N-1)/2, (N-1)/2+1, dtype= int), np.arange(-(N-1)/2, (N-1)/2+1, dtype= int), np.arange(-(N-1)/2, (N-1)/2+1, dtype= int))
	rho_hat = ((-1)**np.abs(wx+wy+wz)/(N**3)*rho_hat)

	#print('Now running Linear Interpolation ...')
	N_range = np.arange(-(N-1)/2, (N-1)/2+1)
	rho_hat = RGI((N_range,N_range,N_range), rho_hat, bounds_error= False, fill_value=0)

	### make 2D mesh that takes slice of rho_hat which NxNx3 grid 
	#print('Creating 2D Sampling grid ...')
	eta_x, eta_y = np.meshgrid(np.arange(-(N-1)/2, (N-1)/2+1, dtype= int), np.arange(-(N-1)/2, (N-1)/2+1, dtype= int), indexing= "ij")
	eta_x= eta_x[...,np.newaxis]
	eta_y= eta_y[...,np.newaxis]

	R= np.transpose(R)
	a= R[0]
	b= R[1]
	grid = eta_x*a + eta_y*b
	sample_grid = rho_hat(grid)

	em_slice = (-1)**np.abs(eta_x[...,0]+eta_y[...,0])*(N**2)*sample_grid

	em_slice = np.fft.ifftn(np.fft.ifftshift(em_slice))

	return np.real(em_slice)

def get_rotation_matrix():
    a = 2*np.random.rand(3)-1
    a /= np.linalg.norm(a)
    b = np.cross(a, 2*np.random.rand(3)-1)
    b /= np.linalg.norm(b)
    c = np.cross(a, b)
    c /= np.linalg.norm(c)

    return np.transpose(np.array([a, b, c]))