#! /usr/bin/env python3
import mrcfile
from scipy.interpolate import RegularGridInterpolator as RGI
import sys
from matplotlib import pyplot as plt
from RotationMatrix import get_rotation_matrix
from Project_fst import project_fst

##########################
#####input
input_file=mrcfile.open('zika_153.mrc')
rotaion_num = int(3)
#####
print ('Now generating random rotation matrix')
Rs = [get_rotation_matrix() for i in range(rotaion_num)]

print ('Now running project_fst')
rho = input_file.data
n = 0
for R in Rs:
	image = project_fst(rho, R)
	plt.imshow(image)
	plt.savefig(str(n)+'.png')
	n += 1
sys.exit()