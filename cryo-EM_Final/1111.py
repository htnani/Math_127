#! /usr/bin/env python3

###
import image_sim as SIM
import reconstruct as RC
import mrcfile
#import common_lines as CL
import os
from PIL import Image
import numpy as np
####
imagesDir = '/Users/Ethan/Desktop/Math_127/cryo-EM/images_69.old'

##This part is uesd for fining the common lines of the two images
print ('Now finding common lines ...') 

im0 = Image.open(imagesDir+'/'+'0.tiff')
im1array = np.array(im)
im1 = Image.open(imagesDir+'/'+'1.tiff')
im1array = np.array(im)
im2 = Image.open(imagesDir+'/'+'2.tiff')
im2array = np.array(im)

def getLine(image, theta, num_points = 100):
	x = np.arange(image.shape[1])
	y = np.arange(image.shape[0])

	if theta == 0:  
		x0 = 0
		y0 = (N - 1)/2  
		x1 = N 
		y1 = (N - 1)/2  
	elif (theta <= np.pi/4 and theta >= -np.pi/4 ) or (theta >= 3*np.pi/4 and theta <= 5*np.pi/4 ):
	 	x0 = 0
	 	y0 = np.tan(theta) * (N - 1)/2 + ((N - 1)/2)
	 	x1 = N
	 	y1 = (N-1)/2 - np.tan(theta) * (N-1)/2
	else:
	 	x0 = np.tan(theta - np.pi/2 ) * (N - 1)/2 + ((N - 1)/2) 
	 	y0 = N
	 	x1 = (N-1)/2 - np.tan(theta - np.pi/2) * (N-1)/2
	 	y1 = 0
	x_values = np.linspace(x0, x1, num_points)
	y_values = np.linspace(y0, y1, num_points)

	z_values = scipy.ndimage.map_coordinates(image, np.vstack((x_values,y_values)))
	# fig, axes = plt.subplots(nrows=2)
	# axes[0].imshow(image)
	# axes[0].plot([x0, x1], [y0, y1], 'ro-')
	# axes[0].axis('image')
	# axes[1].plot(z_values)
	# plt.show()
	return z_values

def commonLine(image1, image2): 
	commonLine_array = []
	commonLine_product_array = []
	for theta_i in range(360):
		zi = getLine(image1, (theta_i * np.pi *2)/360)
		for theta_j in range(360):
			zj = getLine(image2, (theta_j * np.pi *2)/360)
			commonLine_array.append(zi)
			commonLine_array.append(zj)
			commonLine_product_array.append(np.dot(zi,zj))
	print("Done")
	max_index, max_value = max(enumerate(commonLine_product_array), key=lambda p: p[1])
	return commonLine_array[max_index*2], max_value

def getAllLines(image1, image2, image3): 
	line_12 = commonLine(image1, image2)
	line_13 = commonLine(image1, image3)
	line_21 = commonLine(image2, image1)
	line_23 = commonLine(image2, image3)
	line_31 = commonLine(image3, image1)
	line_32 = commonLine(image3, image2)

	return line_12/norm(line_12), line_21/norm(line_21), line_13/norm(line_13),line_23/norm(line_23), line_31/norm(line_31), line_32/norm(line_32)

def getArcLengths(cos_theta):
	return np.arccos(cos_theta/4)

def getGlobalAngles(A, B, C): 
	return np.arccos((np.cos(C)- np.cos(A)*np.cos(B))/(np.sin(A) * np.sin(B)))

#I literally have no idea what to do here 
def getOrientations(some_3D_matrix, angle): 
	a1 = np.array([1,0,0])
	b1 = np.array([0,1,0])
	c1 = np.cross(a1,b1)
	R1 = np.transpose(np.array([a1,b1,c1]))
	
	N = some_3D_matrix.shape[0]
	rot_matrix = np.array([[np.cos(theta), -np.sin(theta), 0], 
							[np.sin(theta), np.cos(theta), 0], 
							[0, 0 , 1]])

	#I think this is the initial alignment step 
	proposed_matrix = np.zeros((N,N,N))
	for theta in range(0,360): 
		theta = theta * np.pi *2/360
		rot_matrix = np.array([[np.cos(theta), -np.sin(theta), 0], 
								[np.sin(theta), np.cos(theta), 0], 
								[0, 0 , 1]])
		proposed_matrix = some_3D_matrix*numpy.linalg.inv(R1)
		if np.isclose(rot_matrix, some_3D_matrix):
			break
	# I think this is the rotation by the angle of interest
	orientation = proposed_matrix * rot_matrix(angle)

	return the_rest_of_the_thing
'''
l12, l21, l13, l23, l31, l32= getAllLines(image1, image2, image3)
C, A, B = getArcLengths(np.dot(l12, l13)), getArcLengths(np.dot(l13, l23)), getArcLengths(np.dot(l12, l23))
gamma, beta, alpha = getGlobalAngles(A, B, C), getGlobalAngles(A, C, B), getGlobalAngles(B, C, A)
'''