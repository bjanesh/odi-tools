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


try:
    object_str, filters, instrument, images, new_extension, remove_tpv_flag, trim_image_flag, wcs_flag, trim_section, airmasses = odi.photcfgparse('phot_confing.yaml')
except IOError:
    print 'phot_config.yaml does not exist, quitting...'
    exit()

images_ = []
new_images_=[]
for filt in images:
    for key in images[filt]:
        images_.append(images[filt][key])
        new_images_.append(images[filt][key].strip('.fits') + new_extension)

nopv_images = new_images_
for i,img in enumerate(images_):
    if trim_image_flag == True:
        x1, x2, y1, y2 = trim_section[0], trim_section[1], trim_section[2], trim_section[3]
        odi.trim_img(img,x1,x2,y1,y2)
    if remove_tpv_flag == True:
        img_nopv = odi.tpv_remove(img.nofits()+'.trim.fits')
        nopv_images[i] = img_nopv

if wcs_flag == True:
    if not os.path.isfile('full_wcs_fix.done'):
        odi.full_sdssmatch(nopv_images[0],nopv_images[1],instrument,gmaglim=23.0)
        odi.fix_wcs_full(nopv_images[0],coords=nopv_images[0][:-5]+'.wcs.coo')
        odi.fix_wcs_full(nopv_images[1],coords=nopv_images[0][:-5]+'.wcs.coo')
        with open('full_wcs_fix.done','w+') as f:
            print >> f, ''

apcor_values = {}
apcor_stds = {}
apcor_sems = {}
fwhm_values = {}
for i,img in enumerate(images_):
    median_fwhm,median_bg_mean,median_bg_median,median_bg_std = odi.read_proc('derived_props.txt',filters[i])
    img = nopv_images[i]
    peaks,gfwhms = odi.getfwhm_full_sdss(img)

    median_gfwhm = np.median(gfwhms[np.where(gfwhms < 20.0)])
    print median_gfwhm
    fwhm_values[i] = median_gfwhm

    odi.sdss_phot_full(img,median_gfwhm,airmasses[i])
    apcor, apcor_std, apcor_sem = odi.apcor_sdss(img, median_gfwhm, inspect=False)
    apcor_values[i] = apcor
    apcor_stds[i] = apcor_std
    apcor_sems[i] = apcor_sem

    odi.find_sources_full(img,median_gfwhm,median_bg_std,threshold=3.5)
    odi.phot_sources_full(img,median_gfwhm,airmasses[i],1.0)
    odi.phot_sources_xy2sky(img,instrument)

photcalFile = odi.calibrate_match(nopv_images[0],nopv_images[1],fwhm_values[0],fwhm_values[1],airmasses[0],airmasses[1])

odi.match_phot_srcs(nopv_images[0],nopv_images[1])
odi.calc_calibrated_mags(apcor_values[0], 0, apcor_values[1], 0, photcalFile, object_str)
