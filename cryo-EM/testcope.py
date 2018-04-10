#! /usr/bin/env python3
import mrcfile
from scipy.interpolate import RegularGridInterpolator as RGI
import sys
from matplotlib import pyplot as plt
import numpy as np
import image_sim as SIM
import reconstruct as RC
import multiprocessing
from functools import partial

##########################
#####input
#input_file=mrcfile.open('emd_8116.map')
input_file=mrcfile.open('zika_153.mrc')
rotaion_num = int(3)

#####
print ('Now generating random rotation matrix ...')
Rs = [SIM.get_rotation_matrix() for i in range(rotaion_num)]
print ('Now running project_fst ...')
rho = input_file.data
n = 0
sim_images_and_Rs = [] #(image, R)

print ('Now producing the images ...')
#multi processing
cores = multiprocessing.cpu_count()
pool = multiprocessing.Pool(processes=cores)

fit_rho = partial(SIM.project_fst, mol = rho)

for image in pool.imap(fit_rho, Rs):
	plt.imshow(image)
	plt.savefig(str(n)+'.png')
	sim_images_and_Rs.append((image, Rs[n]))
	n += 1

print ('Now runing the reconstruct code ...')
N = sim_images_and_Rs[0][0].shape[0] 
B = np.zeros((N,N,N))
L = np.zeros((N,N,N))

for mp_result in pool.imap(RC.reconstruct_mp, sim_images_and_Rs):
	return_B, return_L = mp_result
	B += return_B
	L += return_L

back_img = np.float32(np.real(B)/np.real(L))
back_img = np.nan_to_num(back_img)
print('Now saving the mrc files ...')
with mrcfile.new('full2.mrc', overwrite=True) as mrc:
	mrc.set_data(back_img)

sys.exit()

'''mp_return_list = pool.map(RC.reconstruct_mp, sim_images_and_Rs)
for mp_result in mp_return_list:
	return_B, return_L = mp_result
	B += return_B
	L += return_L'''

#pool.close()

'''back_img = np.float32(np.real(B)/np.real(L))
for 0 in to nan
back_img = np.nan_to_num(back_img)
print('Now saving the mrc files ...')
with mrcfile.new('full2.mrc', overwrite=True) as mrc:
	mrc.set_data(back_img)

sys.exit()'''
