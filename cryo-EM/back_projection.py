from scipy.interpolate import RegularGridInterpolator as RGI
import numpy as np 
import mrcfile

def back_project(R, image):
	N = len(image)
	I_hat = np.fft.fftshift(np.fft.fftn(image))
	rect = np.zeros(N)
	rect[(N-1)/2:N-((N-1)/2)] = 1
	rect_hat = np.fft.fftshift(np.fft.fftn(rect))
	rect_hat = np.tile(rect_hat[np.newaxis, np.newaxis, :], (N, N, 1))
	resutlt = np.zeros((N, N, N))
	I_hat = np.tile(I_hat[..., np.newaxis], (1, 1, N))
	image = np.real(np.fft.ifftn(np.fft.ifftshift(I_hat*rect_hat)))

	####todo: setup density
'''
	ouput = mrcfile.new('output.mrc')
	ouput.set_data(resutlt)'''
