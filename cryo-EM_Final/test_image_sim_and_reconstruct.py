#! /usr/bin/env python3

###
import image_sim as SIM
import reconstruct as RC
import mrcfile
import os
import numpy as np
from PIL import Image

input_file = mrcfile.open('zika_153.mrc')
rotaion_num = int(20)   ##change to how many pictures you need
imagesDir = '/Users/Ethan/Desktop/cryo-EM/images'

#####
##This part is used for generating the images
print ('Now generating random rotation matrix ...')
Rs = [SIM.getRotationMatrix() for i in range(rotaion_num)]

print ('Now running project_fst ...')
rho = input_file.data

##Checking whether the imagesDir exists and saving the simulated images to imagesDir
SIM.checkImagesDir(imagesDir)
images = []
n = 0
#print(Rs)
for R in Rs:
	image = SIM.project_fst(rho, R)
	SIM.saveImage(imagesDir, image, n)
	n += 1
	images.append(image)
	#print(image)

##This part is used for compressing the images with the original orientations 
data = zip(images, Rs)

##This part is used for reconstructing of the 3D model based on the input images and their corresponding orientations
print ('Now running backProjection and reconstruct ...')
sum_b_j_hat, sum_l_j_hat = 0, 0
for image, R in data:
	b_j_hat, l_j_hat = RC.reconstruction(image, R)
	sum_b_j_hat += b_j_hat
	sum_l_j_hat += l_j_hat
back_image = np.fft.ifftn(np.fft.ifftshift(sum_b_j_hat/sum_l_j_hat)).real
#back_image.real

##This part is used for generating the final outputs of the 3D mrc files
print ('Now creating output.mrc ...')
with mrcfile.new('output.mrc', overwrite=True) as mrc:
	mrc.set_data(back_image.astype('float32'))
