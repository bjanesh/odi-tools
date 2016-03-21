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

images_g = glob.glob('*_odi_g*.fits')
images_g.sort()
#print images_g
images_r = glob.glob('*_odi_r*.fits')
images_r.sort()
filters = ['odi_g','odi_r']

images = images_g+images_r


source = 'sdss'
inst = odi.instument(images[0])


g_img = 'GCPair-F6_odi_g.fits'
g_imgr = 'GCPair-F6_odi_g.trim.fits'
odi.trim_img(g_img)
r_img = 'GCPair-F6_odi_r.fits'
r_imgr = 'GCPair-F6_odi_r.trim.fits'
odi.trim_img(r_img)


#Color eq steps
odi.full_sdssmatch(g_img,r_img,inst)

median_fwhm,median_bg_mean,median_bg_median,median_bg_std = odi.read_proc('derived_props.txt','odi_g')
print median_fwhm,median_bg_mean,median_bg_median,median_bg_std
airmass_g = odi.get_airmass(images_g)
odi.sdss_phot_full(g_img,median_fwhm,airmass_g)

median_fwhmr,median_bg_meanr,median_bg_medianr,median_bg_stdr = odi.read_proc('derived_props.txt','odi_r')
print median_fwhmr,median_bg_meanr,median_bg_medianr,median_bg_stdr

airmass_r = odi.get_airmass(images_r)
odi.sdss_phot_full(r_img,median_fwhmr,airmass_r)

odi.calibrate_match(g_img,r_img,median_fwhm,median_fwhmr,airmass_g,airmass_r)

#Phot Steps g
odi.find_sources_full(g_imgr,median_fwhm,median_bg_std)
odi.phot_sources_full(g_imgr,median_fwhm,airmass_g)
odi.phot_sources_xy2sky(g_imgr,inst)

#Phot Steps r
odi.find_sources_full(r_imgr,median_fwhmr,median_bg_stdr)
odi.phot_sources_full(r_imgr,median_fwhmr,airmass_r)
odi.phot_sources_xy2sky(r_imgr,inst)

odi.match_phot_srcs(g_imgr,r_imgr)