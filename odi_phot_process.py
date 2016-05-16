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

# These are not acutally needed.
# We only really need one to get the inst variable

images_g = glob.glob('20*_odi_g*.fits')
images_g.sort()
images_r = glob.glob('20*_odi_r*.fits')
images_r.sort()
filters = ['odi_g','odi_r']

images = images_g+images_r

# Get the airmass of reference image
grefimg = odi.find_ref_image(images_g)
ghduref = odi.fits.open(images_g[grefimg])
airmass_g = ghduref[0].header['airmass']
#airmass_g = 1.045
#airmass_r = 1.161
ghduref.close()

rrefimg = odi.find_ref_image(images_r)
rhduref = odi.fits.open(images_r[rrefimg])
airmass_r = rhduref[0].header['airmass']
rhduref.close()


source = 'sdss'
inst = odi.instrument(images[0])
fitsref = odi.fits.open(images[0])
hduref = fitsref[0]
objname = hduref.header['object']
fitsref.close()
#filter_name = hduref.header['filter']
#sky_med = hduref.header['skybg']
#output = objname+'_'+filter_name+'.fits'

g_img = g_imgr = objname+'_odi_g.fits'
print 'The g image is: ', g_img

r_img = r_imgr = objname+'_odi_r.fits'
print 'The r image is: ', r_img


#Color eq steps
#First find the sdss sources that are in each stacked image
if not os.path.isfile('full_wcs_fix.done'):
    odi.full_sdssmatch(g_img,r_img,inst,gmaglim=19.0)
    odi.fix_wcs_full(g_img,coords=g_img[:-5]+'.wcs.coo')
    odi.fix_wcs_full(r_img,coords=r_img[:-5]+'.wcs.coo')
    with open('full_wcs_fix.done','w+') as f:
        print >> f, ''

odi.full_sdssmatch(g_img,r_img,inst,gmaglim=23.0)
#Get source and background characteristics from 'derived_props.txt'
median_fwhm,median_bg_mean,median_bg_median,median_bg_std = odi.read_proc('derived_props.txt','odi_g')

# Measure gfwhm of sdss stars on combined image
g_peaks,g_gfwhms = odi.getfwhm_full_sdss(g_img)
median_fwhmg = np.median(g_gfwhms[np.where(g_gfwhms < 20.0)])

#Phot sdss g sources on stacked image
odi.sdss_phot_full(g_img,median_fwhmg,airmass_g)


#Repeat same steps above, but for the r stacked image
median_fwhmr,median_bg_meanr,median_bg_medianr,median_bg_stdr = odi.read_proc('derived_props.txt','odi_r')

# Measure gfwhm of sdss stars on combined image
r_peaks,r_gfwhms = odi.getfwhm_full_sdss(r_img)
median_fwhmr = np.median(r_gfwhms[np.where(r_gfwhms < 20.0)])

odi.sdss_phot_full(r_img,median_fwhmr,airmass_r)

#Solve the color equations using the pair of stacked images. This was just
#taken odi_calibrate.
odi.calibrate_match(g_img,r_img,median_fwhmg,median_fwhmr,airmass_g,airmass_r)

apcor_g, apcor_std_g, apcor_sem_g = odi.apcor_sdss(g_img, median_fwhmg,inspect=True)
apcor_r, apcor_std_r, apcor_sem_r = odi.apcor_sdss(r_img, median_fwhmr,inspect=True)

#Phot Steps g
#Find all sources using daofind
odi.find_sources_full(g_imgr,median_fwhmg,median_bg_std,threshold=3.5)
#Phot found sources
odi.phot_sources_full(g_imgr,median_fwhmg,airmass_g,1.0)
#Convert xy positions of sources to Ra Dec
odi.phot_sources_xy2sky(g_imgr,inst)

#Phot Steps r
odi.find_sources_full(r_imgr,median_fwhmr,median_bg_stdr,threshold=3.5)
odi.phot_sources_full(r_imgr,median_fwhmr,airmass_r,1.0)
odi.phot_sources_xy2sky(r_imgr,inst)

#Create a matched catalog of sources on the g and r frame
odi.match_phot_srcs(g_imgr,r_imgr)
odi.calc_calibrated_mags(apcor_g, 0, apcor_r, 0)

