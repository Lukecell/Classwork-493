# -*- coding: utf-8 -*-
"""
Created on Sat May  6 00:01:05 2023

@author: phyfr
"""

import os
from pathlib import Path
from astropy.nddata import CCDData
from astropy.stats import mad_std
import ccdproc as ccdp
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u

data_path = Path(os.getcwd() + '/Documents\ARCSAT_OBSERVATIONS')
calibrated_data = Path(os.getcwd()+ '/Documents\ARCSAT_OBSERVATIONS/Calibrated')
calibrated_data.mkdir(exist_ok=True)


files = ccdp.ImageFileCollection(data_path)

science_images = files.files_filtered(imagetyp = 'Light Frame', include_path = True)
science_size = np.shape(science_images)
# overscan = '[:2048, :]' #PLACEHOLDER OVERSCAN


raw_darks = files.files_filtered(imagetyp = 'Dark Frame', include_path = True)
dark_size = np.shape(raw_darks)
calibrated_darks = raw_darks
for i in range(dark_size[0]):
    darker = CCDData.read(raw_darks[i], unit = 'adu') # PLACEHOLDER UNIT, MIGHT BE ELECTRON OR PHOTON
    # overscan_subtracted = ccdp.subtract_overscan(darker, overscan = darker[:, :2048])
    # trimmer = ccdp.trim_image(overscan_subtracted[:, :2048])
   
    calibrated_darks[i] = darker #trimmer if overscan exists

combined_darks = ccdp.combine(calibrated_darks,
                              method='average',
                              sigma_clip=True, sigma_clip_low_thresh=3, sigma_clip_high_thresh=5,
                              sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std,
                              mem_limit=350e8
                            )

combined_darks.meta['combined'] = True

combined_darks.write(calibrated_data / 'combined_darks.fit', overwrite = True)

###########################################################################################################################################################################

raw_flats = files.files_filtered(imagetyp = 'FLAT', include_path = True)

filters = ['sdss_u','sdss_g','sdss_r','sdss_i','sdss_z',]



variables = {filters[i]: files.files_filtered(imagetyp = 'FLAT',FILTER = filters[i], include_path = True) for i in range(len(filters))}
for i in range(len(variables)):
    if np.size(variables[filters[i]]) == 0:
        del variables[filters[i]]


    

combed_filters = []
for i in range(len(variables)):
    combed_filters.append(list(variables.keys())[i])
    
flat_u = []
flat_g = []
flat_r = []
flat_i = []
flat_z = []
    
for i in range(len(combed_filters)):
    
    if combed_filters[i] == 'sdss_u':
        flat_u = list(variables.values())[i]
        
        for j in range(len(flat_u)):
            flat_u[j] = CCDData.read(flat_u[j], unit = 'adu')
            
    if combed_filters[i] == 'sdss_g':
        flat_g = list(variables.values())[i]

        for j in range(len(flat_g)):
            flat_g[j] = CCDData.read(flat_g[j], unit = 'adu')
            
    if combed_filters[i] == 'sdss_r':
        flat_r = list(variables.values())[i]

        for j in range(len(flat_r)):
            flat_r[j] = CCDData.read(flat_r[j], unit = 'adu')
            
    if combed_filters[i] == 'sdss_i':
        flat_i = list(variables.values())[i]

        for j in range(len(flat_i)):
            flat_i[j] = CCDData.read(flat_i[j], unit = 'adu')
            
    if combed_filters[i] == 'sdss_z':
        flat_g = list(variables.values())[i]

        for j in range(len(flat_z)):
            flat_z[j] = CCDData.read(flat_z[j], unit = 'adu')
        

combined_flats = []

for i in range(len(combed_filters)):
    
    combined_flatser =  ccdp.combine(variables[combed_filters[i]],
                                  method='average',
                                  sigma_clip=True, sigma_clip_low_thresh=3, sigma_clip_high_thresh=5,
                                  sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std,
                                  mem_limit=350e8
                                )
    combined_flats.append(combined_flatser)
    combined_flats[i].meta['combined'] = True

    combined_flats[i].write(os.getcwd() + '/Documents\ARCSAT_OBSERVATIONS/Calibrated/combined_flats'+ combed_filters[i] +'.fit', overwrite = True)



u_science = files.files_filtered(imagetyp = 'Light Frame', FILTER = 'sdss_u', include_path = True)
g_science = files.files_filtered(imagetyp = 'Light Frame', FILTER = 'sdss_g', include_path = True)
r_science = files.files_filtered(imagetyp = 'Light Frame', FILTER = 'sdss_r', include_path = True)
i_science = files.files_filtered(imagetyp = 'Light Frame', FILTER = 'sdss_i', include_path = True)
z_science = files.files_filtered(imagetyp = 'Light Frame', FILTER = 'sdss_z', include_path = True)

for i in range(len(u_science)):
    u_science[i] = CCDData.read(u_science[i], unit = 'adu')

for i in range(len(g_science)):
    g_science[i] = CCDData.read(g_science[i], unit = 'adu')
    
for i in range(len(r_science)):
    r_science[i] = CCDData.read(r_science[i], unit = 'adu')

for i in range(len(i_science)):
    i_science[i] = CCDData.read(i_science[i], unit = 'adu')

for i in range(len(z_science)):
    z_science[i] = CCDData.read(z_science[i], unit = 'adu')


m13_g = g_science[0]

m13_g = ccdp.ccd_process(m13_g, dark_frame = combined_darks, dark_exposure = 1800*u.s, data_exposure = 120*u.s, master_flat = combined_flats[0])
m13_g.write(os.getcwd() + '/Documents\ARCSAT_OBSERVATIONS/Calibrated/M13_g.fit', overwrite = True)

NGC_4726_g = g_science[1]

NGC_4726_g = ccdp.ccd_process(NGC_4726_g, dark_frame = combined_darks, dark_exposure = 1800*u.s, data_exposure = 900*u.s, master_flat = combined_flats[0])
NGC_4726_g.write(os.getcwd() + '/Documents\ARCSAT_OBSERVATIONS/Calibrated/NGC_4726_g.fit', overwrite = True)

plt.imshow(m13_g)

# combed_science = []
# final_images = science_images
# read_science = science_images
# for i in range(science_size[0]):
#     read_science[i] = CCDData.read(science_images[i], unit = 'adu')
# plt.imshow(read_science[1])
# plt.show()
# for i in range(science_size[0]):
    
#     final_images[i] = ccdp.ccd_process(read_science[i], dark_frame = combined_darks, dark_exposure = 1800*u.s, data_exposure = 300*u.s, master_flat = combined_flats)
#     final_images[i].write(os.getcwd() + '/Downloads/fitsers/example-cryo-LFC/calibrated/' + + str([i]) +'.fit', overwrite = True)
    
# plt.imshow(final_images[1])
# NEEDED IF DARKS DONT INCLUDE BIAS 
# raw_bias = files.files_filtered(imagetyp = 'BIAS', include_path = True)
# size = np.shape(raw_bias)
# calibrated_biases = raw_bias
# for i in range(size[0]):
#     biaser = CCDData.read(raw_bias[i], unit = 'adu') # PLACEHOLDER UNIT, MIGHT BE ELECTRON OR PHOTON
#     overscan_subtracted = ccdp.subtract_overscan(biaser, overscan = biaser[:, :2048])
#     trimmer = ccdp.trim_image(overscan_subtracted[:, :2048])
    
#     calibrated_biases[i] = trimmer

# combined_bias = ccdp.combine(calibrated_biases,
#                               method='average',
#                               sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
#                               sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std,
#                               mem_limit=350e6
#                             )

# combined_bias.meta['combined'] = True

# combined_bias.write(calibrated_data / 'combined_bias.fit', overwrite = True)

