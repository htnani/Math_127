#! /usr/bin/env python3

###
import image_sim as SIM
import reconstruct as RC
import mrcfile
import common_lines as CL
from PIL import Image
####
import sys
import numpy as np
import reconstruct as RC
import multiprocessing
from functools import partial

input_file=mrcfile.open('zika_153.mrc')
rotaion_num = int(10)

#####
print ('Now generating random rotation matrix ...')
Rs = [SIM.getRotationMatrix() for i in range(rotaion_num)]
print (Rs)
print ('Now running project_fst ...')
rho = input_file.data

images = []
'''for R in Rs:
	image = SIM.project_fst(rho, R)
	SIM.saveImage(image, n)
	n += 1
	images.append(image)

data = zip(images, Rs)

print ('Now running backProjection and reconstruct ...')
sum_b_j_hat, sum_l_j_hat = 0, 0
for image, R in data:
	b_j_hat, l_j_hat = RC.reconstruction(image, R)
	sum_b_j_hat += b_j_hat
	sum_l_j_hat += l_j_hat
back_image = np.fft.ifftn(np.fft.ifftshift(sum_b_j_hat/sum_l_j_hat)).real
#back_image.real


print ('Now creating full2.mrc ...')
with mrcfile.new('full2.mrc', overwrite=True) as mrc:
	mrc.set_data(back_image.astype('float32'))
'''
print ('Now finding common lines ...') 

im = Image.open('a_image.tif')
im.show()