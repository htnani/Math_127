import numpy as np
import mrcfile
from scipy.interpolate import RegularGridInterpolator as RGI

def reconstruct_mp(image_array_and_R_array): #(image_array, R_array)
	image_array, R_array = image_array_and_R_array
	N = image_array.shape[0] 

	#print("Computing image_hat")
	image_hat = np.fft.fftn(image_array)
	image_hat = np.fft.fftshift(image_hat)
	freq_range = np.arange(-(N-1)/2, (N-1)/2 +1, dtype= int) 
	wx_j, wy_j = np.meshgrid(freq_range, freq_range)
	image_hat = (((-1)**np.abs(wx_j+wy_j))/N**2)*image_hat
	image_hat = np.tile(image_hat[...,np.newaxis], [1,1,N])

	#print("Computing l_j_hat")
	wz_j = freq_range
	l_j = N * np.sinc(N * np.pi * wz_j)
	l_j_hat = np.tile(l_j[np.newaxis,np.newaxis,...],[N,N,1])

	#print("Creating rotated grid")
	wlx, wly, wlz = np.meshgrid(freq_range, freq_range, freq_range)

	rot_grid = np.zeros((N,N,N,3))

	#todo: use matrix?
	for j in range(N):
		for k in range(N):
			for l in range(N):
				point = np.array([wlx[j,k,l], wly[j,k,l], wlz[j,k,l]]).reshape(3,1)
				rot_grid[j,k,l] = np.dot(np.transpose(R_array), point).reshape(1,3) ###why this is diffrent from the note? i chenged this

	image_hat_fnc = RGI((freq_range,freq_range,freq_range), image_hat, bounds_error= False, fill_value=0)
	b_j_hat_fnc = RGI((freq_range,freq_range,freq_range), l_j_hat, bounds_error= False, fill_value=0)

	image_hat = image_hat_fnc(rot_grid)
	l_j_hat = image_hat_fnc(rot_grid)

	print (l_j_hat)

	b_j_hat = image_hat * l_j_hat
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