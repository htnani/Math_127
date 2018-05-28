from scipy.interpolate import RegularGridInterpolator as RGI
import numpy as np 
from PIL import Image
import os, random


def getRotationMatrix():
	'''
	this function return a random normalized 3D vector 
	'''
	a = 2*np.random.rand(3)-1
	a /= np.linalg.norm(a)
	b = np.cross(a, 2*np.random.rand(3)-1)
	b /= np.linalg.norm(b)
	c = np.cross(a, b)
	return np.transpose(np.array([a, b, c]))


def project_fst(mol, R):
	'''
	mol is an NXNXN array that contains the electron density of the protein thing
	R is the rotation matrix that symbolizes the viewing direction of the microscope
	this function is used to inicialize the molecular and return the em slice of the molecular
	'''
	
	#print('Now running 3D fourier Transform ...')
	N = len(mol)
	NRange = np.arange(-(N-1)/2, (N-1)/2+1, dtype = int)
	rhoHat = np.fft.fftshift(np.fft.fftn(mol))
	wx, wy, wz = np.meshgrid(NRange, NRange, NRange)
	rhoHat = ((-1)**np.abs(wx+wy+wz)/(N**3)*rhoHat)
	#print('Now running Linear Interpolation ...')
	rhoHat = RGI((NRange, NRange, NRange), rhoHat, bounds_error= False, fill_value=0)

	### make 2D mesh that takes slice of rhoHat which NxNx3 grid 
	#print('Creating 2D Sampling grid ...')
	etaX, etaY = np.meshgrid(NRange, NRange, indexing= "ij")
	etaX= etaX[...,np.newaxis]
	etaY= etaY[...,np.newaxis]

	R= np.transpose(R)
	a= R[0]
	b= R[1]
	grid = etaX*a + etaY*b
	sampleGrid = rhoHat(grid)

	emSlice = (-1)**np.abs(etaX[...,0]+etaY[...,0])*(N**2)*sampleGrid

	emSlice = np.fft.ifftn(np.fft.ifftshift(emSlice))

	return np.real(emSlice)


def saveImage(imagesDir, image, seq):
	#image = SIM.project_fst(Rs[0], rho)
	im = Image.fromarray(image)
	im.save(imagesDir + '/' + str(seq) + '.tiff')


def checkImagesDir(imagesDir):
	if not os.path.exists(imagesDir):
		os.makedirs(imagesDir)
	else:
		os.system('mv '+imagesDir+' '+imagesDir+'_'+str(random.randint(0, 100))+'.old')
		os.makedirs(imagesDir)
	print('Saving the simulated images to ', imagesDir)
