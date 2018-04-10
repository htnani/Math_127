import numpy as np
import mrcfile
from scipy.interpolate import RegularGridInterpolator as RGI

##
from scipy import fftpack, interpolate, linalg, ndimage

def reconstructdie(imageArray_and_RArray): ##back projections
	#(imageArray, RArray)
	imageArray, RArray = imageArray_and_RArray
	N = imageArray.shape[0] 

	#print('Computing imageHat')
	imageHat = np.fft.fftshift(np.fft.fftn(imageArray))
	NRange = np.arange(-(N-1)/2, (N-1)/2+1, dtype= int) 
	wxJ, wyJ = np.meshgrid(NRange, NRange)
	imageHat = (((-1)**np.abs(wxJ+wyJ))/N**2)*imageHat
	imageHat = np.tile(imageHat[...,np.newaxis], [1,1,N])

	#print("Computing l_j_hat")
	wz_j = NRange
	l_j = N * np.sinc(N * np.pi * wz_j)
	l_j_hat = np.tile(l_j[np.newaxis,np.newaxis,...],[N,N,1])

	#print("Creating rotated grid")
	wlx, wly, wlz = np.meshgrid(NRange, NRange, NRange)

	rot_grid = np.zeros((N,N,N,3))

	#todo: use matrix?
	for j in range(N):
		for k in range(N):
			for l in range(N):
				point = np.array([wlx[j,k,l], wly[j,k,l], wlz[j,k,l]]).reshape(3,1)
				rot_grid[j,k,l] = np.dot(np.transpose(RArray), point).reshape(1,3) ###why this is diffrent from the note? i chenged this

	imageHat_fnc = RGI((NRange,NRange,NRange), imageHat, bounds_error= False, fill_value=0)
	b_j_hat_fnc = RGI((NRange,NRange,NRange), l_j_hat, bounds_error= False, fill_value=0)

	imageHat = imageHat_fnc(rot_grid)
	l_j_hat = imageHat_fnc(rot_grid)

	print (l_j_hat)

	b_j_hat = imageHat * l_j_hat
	b_j_hat = ((-1)**np.abs(wlx+wly+wlz)*(N**3))*b_j_hat
	b_j_hat = np.fft.ifftshift(b_j_hat)
	b_j = np.fft.ifftn(b_j_hat)
	
	return (np.real(b_j), np.real(l_j_hat))


#print(back_img.shape)
#with file('back_img.mrc', 'w') as outfile:
#    for slice_2d in image:
#        np.savetxt(outfile, slice_2d)
#np.savetxt('back_img.mrc', back_img, delimiter=',')
#back_obj = open("back_img.mrc", "w")
#back_obj.write(back_img)
#back_obj.close()






##Done  #bp and reconstruct
def backProjection(imageArray, RArray): ##back projections
	N = imageArray.shape[0] 
	#print('Computing imageHat')
	imageHat = np.fft.fftshift(np.fft.fftn(imageArray))
	NRange = np.arange(-(N-1)/2, (N-1)/2+1, dtype= int) 
	wxJ, wyJ = np.meshgrid(NRange, NRange)
	imageHat = (((-1)**np.abs(wxJ+wyJ))/N**2)*imageHat
	imageHat = np.tile(imageHat[...,np.newaxis], [1,1,N])

	#print("Computing l_j_hat")
	wz_j = NRange
	l_j = N * np.sinc(N * np.pi * wz_j)
	l_j_hat = np.tile(l_j[np.newaxis, np.newaxis,...],[N,N,1])

	#print("Creating rotated grid")
	wlx, wly, wlz = np.meshgrid(NRange, NRange, NRange)
	rot_grid = np.zeros((N,N,N,3))   #this is a 3D 0 array

	for j in range(N):
		for k in range(N):
			for l in range(N):
				point = np.array([wlx[j,k,l], wly[j,k,l], wlz[j,k,l]]).reshape(3,1)
				rot_grid[j,k,l] = np.dot(np.transpose(RArray), point).reshape(1,3) 
				'''point = np.array([wlx[j,k,l], wly[j,k,l], wlz[j,k,l]]) 
																rot_grid[j,k,l] = np.dot(RArray, point)
																print (rot_grid[j,k,l])
																print ('ok')'''
	imageHat_fnc = RGI((NRange, NRange, NRange), imageHat, bounds_error= False, fill_value=0)
	b_j_hat_fnc = RGI((NRange, NRange, NRange), l_j_hat, bounds_error= False, fill_value=0)

	imageHat = imageHat_fnc(rot_grid)
	l_j_hat = imageHat_fnc(rot_grid)

	b_j_hat = imageHat * l_j_hat
	b_j_hat = ((-1)**np.abs(wlx+wly+wlz)*(N**3))*b_j_hat
	b_j = np.fft.ifftn(np.fft.ifftshift(b_j_hat))
	
	return np.real(b_j)


##Done
def backProjection2(imageArray): ##back projections
	N = imageArray.shape[0] 
	#print('Computing imageHat')
	imageHat = np.fft.fftshift(np.fft.fftn(imageArray))
	# scale
	NRange = np.arange(-(N-1)/2, (N-1)/2+1, dtype= int)
	etax, etay = np.meshgrid(NRange, NRange, indexing='ij')
	imageHat = np.exp(np.pi * 1j * (etax + etay)) * N**2 * imageHat
	
	# compute F{l_j}
	wz = NRange
	l_j = np.sinc(wz)  ###???
	
	# tiling
	l_j_hat = np.tile(l_j[np.newaxis, np.newaxis,...],[N,N,1])   ###???
	imageHat = np.tile(imageHat[..., np.newaxis], (1, 1, N))
	
	b_j_hat = imageHat * l_j_hat
	return b_j_hat, np.sum(l_j_hat)

##Done
def reconstruction(image, R):
	R= np.transpose(R)
	a= R[0]
	b= R[1]

	b_j_hat, l_j_hat = backProjection2(image)
	
	N = image.shape[0]
	NRange = np.arange(-(N-1)/2, (N-1)/2+1, dtype= int)
	# rotated coordinates
	rot_x, rot_y, rot_z = np.meshgrid(np.arange(-(N-1)/2 , (N-1)/2 + 1, dtype=int), 
							np.arange(-(N-1)/2 , (N-1)/2 + 1, dtype=int),
							np.arange(-(N-1)/2 , (N-1)/2 + 1, dtype=int))
		
	rot_x = rot_x[..., np.newaxis] * a
	rot_y = rot_y[..., np.newaxis] * b
	rot_z = rot_z[..., np.newaxis] * np.cross(a, b)
		
	grid = rot_x + rot_y + rot_z
		
	# standard coordinates
	x = np.arange(-(N-1)/2, (N-1)/2+1)
	y = np.arange(-(N-1)/2, (N-1)/2+1)
	z = np.arange(-(N-1)/2, (N-1)/2+1)
		
	b_j_hat_f = (interpolate.RegularGridInterpolator(points=[x, y, z], 
													 values=b_j_hat, bounds_error=False, fill_value=0))
		
	interpolated = b_j_hat_f(grid)
	w_x, w_y, w_z = np.meshgrid(np.arange(-(N-1)/2 , (N-1)/2 + 1, dtype=int), np.arange(-(N-1)/2 , (N-1)/2 + 1, dtype=int), np.arange(-(N-1)/2 , (N-1)/2 + 1, dtype=int))
	interpolated = np.exp(-np.pi * 1j * (w_x + w_y + w_z)) / N**3 * interpolated
	return interpolated, l_j_hat