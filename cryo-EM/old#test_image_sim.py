#! /usr/bin/env python3
import mrcfile
from scipy.interpolate import RegularGridInterpolator as RGI
import sys
from matplotlib import pyplot as plt
import image_sim as SIM
import multiprocessing
from functools import partial

##########################
#####input
#input_file=mrcfile.open('emd_8116.map')
input_file=mrcfile.open('zika_153.mrc')
rotaion_num = int(20)

#####
print ('Now generating random rotation matrix ...')
Rs = [SIM.get_rotation_matrix() for i in range(rotaion_num)]
print ('Now running project_fst ...')
rho = input_file.data
n = 0
sim_images_list = []

#multi processing
cores = multiprocessing.cpu_count()
pool = multiprocessing.Pool(processes=cores)

fit_rho = partial(SIM.project_fst, mol = rho)

for image in pool.imap(fit_rho, Rs):
	plt.imshow(image)
	plt.savefig(str(n)+'.png')
	n += 1
	sim_images_list.append(image)

'''
#single processing
for R in Rs:
	image = SIM.project_fst(R, rho)
	plt.imshow(image)
	plt.savefig(str(n)+'.png')
	n += 1
	'''

sys.exit()

