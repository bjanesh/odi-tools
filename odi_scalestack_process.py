#!/usr/bin/env python

import sys, os, glob, string
import numpy as np
import astropy as ast
import matplotlib.pyplot as plt
from pyraf import iraf
from tqdm import tqdm
import odi_config as odi
import glob
import shutil
import pandas as pd

images_g = glob.glob('*_KPNO_Ha_8nm*.fits')
images_g.sort()
#print images_g
images_r = glob.glob('*_odi_r*.fits')
images_r.sort()
filters = ['KPNO_Ha_8nm','odi_r']
# 
images = images_g+images_r

# filters = ['odi_r']
# images = images_r

rad, decd = odi.get_targ_ra_dec(images[0], 'OTA33.SCI')

source = 'sdss'
inst = odi.instrument(images[0])

# Scaling with all sources
# for img in images_g:
#     dither  = img.split('.')[1][0]+'_'
#     for key in tqdm(odi.OTA_dictionary):
#         ota = odi.OTA_dictionary[key]
#         if not os.path.isfile(odi.sourcepath+'source_'+ota+'.'+str(img[16:-5])+'.csv'):
#             odi.source_find(img,ota,inst)
#             gaps = odi.get_gaps_rep(img, ota)
#             odi.source_xy(img,ota,gaps,filters[0],inst)
#             fwhm = odi.getfwhm_source(img,ota)
#             print fwhm
#             odi.phot_sources(img, ota, fwhm)
#             odi.phot_combine(img, ota)
#     if not os.path.isfile(odi.sourcepath+dither+filters[0]+'.allsource'):
#         dither_total = odi.sourcepath+dither+filters[0]+'.allsource' 
#         cat_command = 'cat' + ' ' + 'sources/'+'*'+dither+'*_'+filters[0]+'*.totphot' + '>'+dither_total
#         os.system(cat_command)
# 
# refimg_g = odi.find_ref_image(images_g)
# ref_img = images_g[refimg_g]
# 
# scales_g = {}
# stds_g = {}
# n_g = {}
# for img in images_g:
#     scale,std,n = odi.source_scale(img,ref_img,filters[0])
#     scales_g[img] = scale
#     stds_g[img] = std
#     n_g[img] = n
# 
# for img in images_g:
#     for key in odi.OTA_dictionary:
#         ota = odi.OTA_dictionary[key]
#         if not os.path.isfile(odi.scaledpath+'scaled_'+ota+'.'+str(img[16:])):
#             gaps = odi.get_gaps_rep(img, ota)
#             odi.scale_ota(img, ota, scales_g[img])
#             odi.force_update_bpm(img, ota)
# 
# g_img = odi.stack_images(ref_img)

for img in images_r:
    dither  = img.split('.')[1][0]+'_'
    for key in tqdm(odi.OTA_dictionary):
        ota = odi.OTA_dictionary[key]
        if not os.path.isfile(odi.sourcepath+'source_'+ota+'.'+str(img[16:-5])+'.csv'):
            odi.source_find(img,ota,inst)
            gaps = odi.get_gaps_rep(img, ota)
            odi.source_xy(img,ota,gaps,filters[1],inst)
            fwhm = odi.getfwhm_source(img,ota)
            odi.phot_sources(img, ota, fwhm)
            odi.phot_combine(img, ota)
    if not os.path.isfile(odi.sourcepath+dither+filters[1]+'.allsource'):
        dither_total = odi.sourcepath+dither+filters[1]+'.allsource' 
        cat_command = 'cat' + ' ' + 'sources/'+'*'+dither+'*_'+filters[1]+'*.totphot' + '>'+dither_total
        os.system(cat_command)


refimg_r = odi.find_ref_image(images_r)
ref_img = images_r[refimg_r]

scales_r = {}
stds_r = {}
for img in images_r:
    scale,std,n = odi.source_scale(img,ref_img,filters[1])
    scales_r[img] = scale

    
for img in images_r:
    for key in odi.OTA_dictionary:
        ota = odi.OTA_dictionary[key]
        if not os.path.isfile(odi.scaledpath+'scaled_'+ota+'.'+str(img[16:])):
            gaps = odi.get_gaps_rep(img, ota)
            odi.scale_ota(img, ota, scales_r[img])
            odi.force_update_bpm(img, ota)

r_img = odi.stack_images(ref_img)