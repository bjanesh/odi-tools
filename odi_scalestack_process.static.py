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

images_g = glob.glob('*_odi_NB695*.fits')
images_g.sort()
#print images_g
images_r = glob.glob('*_odi_r*.fits')
images_r.sort()
#images_i = glob.glob('*_odi_i*.fits')
#images_i.sort()
filters = ['odi_NB695','odi_r']
# filters = ['odi_r']

images = images_g+images_r

rad, decd = odi.get_targ_ra_dec(images[0], 'OTA33.SCI')

source = 'sdss'
inst = odi.instrument(images[0])

imgnum,fwhm_d,zp_med, zp_std, bg_mean, bg_median, bg_std = np.loadtxt('derived_props.txt',usecols=(0,3,4,5,6,7,8),unpack=True)
ota_d, filt_d = np.loadtxt('derived_props.txt',usecols=(1,2),unpack=True,dtype=str)
id_d = zip(imgnum,ota_d,filt_d)
fwhm_dict = dict(zip(id_d,fwhm_d))

# Scaling with all sources
for img in images_g:
    dither  = img.split('.')[1][0]+'_'
    for key in tqdm(odi.OTA_dictionary):
        ota = odi.OTA_dictionary[key]
        img_id = (int(img.split('.')[1][0]),ota,'odi_NB695')
        if not os.path.isfile(odi.sourcepath+'source_'+ota+'.'+img.base()+'.csv'):
            odi.source_find(img,ota,inst)
            gaps = odi.get_gaps_rep(img, ota)
            odi.source_xy(img,ota,gaps,filters[0],inst)
            fwhm = odi.getfwhm_source(img,ota)
            #fwhm = fwhm_dict[img_id]
            print fwhm
            odi.phot_sources(img, ota, fwhm)
            odi.phot_combine(img, ota)
    if not os.path.isfile(odi.sourcepath+dither+filters[0]+'.allsource'):
        dither_total = odi.sourcepath+dither+filters[0]+'.allsource' 
        cat_command = 'cat' + ' ' + 'sources/'+'*'+dither+'*_'+filters[0]+'*.totphot' + '>'+dither_total
        os.system(cat_command)

refimg_g = odi.find_ref_image(images_g)
ref_img = images_g[refimg_g]

scales_g = {}
stds_g = {}
n_g = {}
for img in images_g:
    scale,std,n = odi.source_scale(img,ref_img,filters[0])
    scales_g[img] = scale
    stds_g[img] = std
    n_g[img] = n

ims = scales_g.keys()
scls = scales_g.values()
new_ref = ims[np.argmax(scls)]
if new_ref != ref_img:  
    ref_img = new_ref  
    for img in images_g:
        scale,std,n = odi.source_scale(img,ref_img,filters[0])
        scales_g[img] = scale
        stds_g[img] = std
        n_g[img] = n
    
for img in images_g:
    for key in odi.OTA_dictionary:
        ota = odi.OTA_dictionary[key]
        if not os.path.isfile(odi.scaledpath+'scaled_'+ota+'.'+img.stem()):
            gaps = odi.get_gaps_rep(img, ota)
            odi.scale_ota(img, ota, scales_g[img])
            odi.force_update_bpm(img, ota)

g_img = odi.stack_images(ref_img)

for img in images_r:
    dither  = img.split('.')[1][0]+'_'
    for key in tqdm(odi.OTA_dictionary):
        ota = odi.OTA_dictionary[key]
        if not os.path.isfile(odi.sourcepath+'source_'+ota+'.'+img.base()+'.csv'):
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

ims = scales_r.keys()
scls = scales_r.values()
new_ref = ims[np.argmax(scls)]
print ref_img, new_ref
if new_ref != ref_img:  
    ref_img = new_ref  
    for img in images_r:
        scale,std,n = odi.source_scale(img,ref_img,filters[0])
        scales_r[img] = scale
        stds_r[img] = std

    
for img in images_r:
    for key in odi.OTA_dictionary:
        ota = odi.OTA_dictionary[key]
        if not os.path.isfile(odi.scaledpath+'scaled_'+ota+'.'+img.stem()):
            gaps = odi.get_gaps_rep(img, ota)
            odi.scale_ota(img, ota, scales_r[img])
            odi.force_update_bpm(img, ota)

r_img = odi.stack_images(ref_img)


for img in images_i:
    dither  = img.split('.')[1][0]+'_'
    for key in tqdm(odi.OTA_dictionary):
        ota = odi.OTA_dictionary[key]
        if not os.path.isfile(odi.sourcepath+'source_'+ota+'.'+img.base()+'.csv'):
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


refimg_i= odi.find_ref_image(images_i)
ref_img = images_i[refimg_i]

scales_i = {}
stds_i = {}
for img in images_i:
    scale,std,n = odi.source_scale(img,ref_img,filters[1])
    scales_i[img] = scale


for img in images_i:
    for key in odi.OTA_dictionary:
        ota = odi.OTA_dictionary[key]
        if not os.path.isfile(odi.scaledpath+'scaled_'+ota+'.'+img.stem()):
            gaps = odi.get_gaps_rep(img, ota)
            odi.scale_ota(img, ota, scales_r[img])
            odi.force_update_bpm(img, ota)

i_img = odi.stack_images(ref_img)
