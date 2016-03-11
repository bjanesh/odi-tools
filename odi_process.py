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

images_g = glob.glob('*GC*_odi_g*.fits')
images_g.sort()
#print images_g
images_r = glob.glob('*GC*_odi_r*')
images_r.sort()
filters = ['odi_g','odi_r']

images = images_g+images_r

rad, decd = odi.get_targ_ra_dec(images[0], 'OTA33.SCI')

# Create offlines catalogs
for img in images:
    print 'Retrieving QR SDSS and 2MASS catalogs for:', img
    for key in odi.OTA_dictionary:
	ota = odi.OTA_dictionary[key]
	outputsd = odi.sdsspath+'offline_'+ota+'.'+str(img[16:-5])+'.sdss'
	x,y = odi.get_sdss_coords_offline(img,ota,output=outputsd)
	output2m = odi.twomasspath+'offline_'+ota+'.'+str(img[16:-5])+'.mass'
	x,y = odi.get_2mass_coords_offline(img,ota,output=output2m)

listfiles = glob.glob('*.lis')
if len(listfiles) == 0:
    odi.imcombine_lists(images, filters)
else:
    print 'imcombine lists done'

for img in images:
    print 'updating bpms for', img
    for key in tqdm(odi.OTA_dictionary):
        ota = odi.OTA_dictionary[key]
        odi.make_bpms(img, ota)

listfiles = glob.glob(odi.skyflatpath+'*.med.fits')
if len(listfiles) == 0:
    for filter in filters:
        odi.dark_sky_flat(filter)
else:
    print 'dark sky flats done'

source = 'sdss'
#source = 'twomass'
if not os.path.isfile('derived_props.txt'):
    f1 = open('derived_props.txt','w+')
    print >> f1, '# img  ota  filter fwhm  zp_med  zp_std  bg_mean  bg_med  bg_std'
    for img in images:
        for key in tqdm(odi.OTA_dictionary):
            ota = odi.OTA_dictionary[key]
            # listfiles = glob.glob(illcorpath+'*.fits')
            # if len(listfiles) == 0:
            image_to_correct = img+'['+ota+']'
            hdulist = odi.fits.open(img)
            hdr = hdulist[0].header
            filt = hdr['filter']
            correction_image = ota+'.'+filt+'.med.fits'
            corrected_image = 'illcor_'+ota+'.'+str(img[16:])
            if not os.path.isfile(odi.illcorpath+corrected_image):
                odi.illumination_corrections(image_to_correct, correction_image, corrected_image)
            # listfiles = glob.glob(reprojpath+'*.fits')
            # if len(listfiles) == 0:
            gaps = odi.get_gaps(img, ota)
            reprojed_image = 'reproj_'+ota+'.'+str(img[16:])
            if not os.path.isfile(odi.reprojpath+reprojed_image):
                pixcrd3 = odi.list_wcs_coords(img, ota, gaps, output=img[:-5]+'.'+ota+'.radec.coo', gmaglim=23., stars_only=True, offline = True, source = source)
                odi.fix_wcs(img, ota, coords=img[:-5]+'.'+ota+'.radec.coo', iters=3)
                odi.reproject_ota(img, ota, rad, decd)
            gaps = odi.get_gaps_rep(img, ota)
            odi.refetch_sdss_coords(img, ota, gaps,gmaglim=20.,offline = True,source=source)
            #run an additional refetch to get the xy for 2mass so they can be used for scaling
            odi.repoxy_offline(img, ota, gaps, gmaglim=20.,source='twomass')
            fwhm = odi.getfwhm_ota(img, ota)
            if source == 'sdss':
		zp_med, zp_std, phot_tbl = odi.zeropoint_ota(img, ota, fwhm)
	    if source == 'twomass':
		zp_med, zp_std = 99.99,99.99
            if not os.path.isfile(odi.bgsubpath+'bgsub_'+ota+'.'+str(img[16:])):
                bg_mean, bg_median, bg_std = odi.bgsub_ota(img, ota, apply=True)
            else:
                bg_mean, bg_median, bg_std = odi.bgsub_ota(img, ota, apply=False)
	    print >> f1, img[16], ota, filt, fwhm, zp_med, zp_std, bg_mean, bg_median, bg_std
    f1.close()


### Scaling with all sources
for img in images_g:
    dither  = img.split('.')[1][0]+'_'
    for key in tqdm(odi.OTA_dictionary):
	ota = odi.OTA_dictionary[key]
	if not os.path.isfile(odi.sourcepath+'source_'+ota+'.'+str(img[16:-5])+'.csv'):
	    odi.source_find(img,ota)
	    gaps = odi.get_gaps_rep(img, ota)
	    odi.source_xy(img,ota,gaps,filters[0])
	    fwhm = odi.getfwhm_source(img,ota)
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

for img in images_g:
    for key in odi.OTA_dictionary:
        ota = odi.OTA_dictionary[key]
        if not os.path.isfile(odi.scaledpath+'scaled_'+ota+'.'+str(img[16:])):
            gaps = odi.get_gaps_rep(img, ota)
            odi.scale_ota(img, ota, scales_g[img])
            odi.force_update_bpm(img, ota)
g_img = odi.stack_images(ref_img)

for img in images_r:
    dither  = img.split('.')[1][0]+'_'
    for key in tqdm(odi.OTA_dictionary):
	ota = odi.OTA_dictionary[key]
	if not os.path.isfile(odi.sourcepath+'source_'+ota+'.'+str(img[16:-5])+'.csv'):
	    odi.source_find(img,ota)
	    gaps = odi.get_gaps_rep(img, ota)
	    odi.source_xy(img,ota,gaps,filters[1])
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
