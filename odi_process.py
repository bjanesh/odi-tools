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
import time

images_g = glob.glob('*_odi_g*.fits')
images_g.sort()
#print images_g
images_r = glob.glob('*_odi_r*.fits')
images_r.sort()
# images_i = glob.glob('*_odi_i*.fits')
# images_i.sort()
filters = ['odi_g','odi_r']#,'odi_i']

images = images_g+images_r#+images_i

rad, decd = odi.get_targ_ra_dec(images[0], 'OTA33.SCI')

source = 'sdss'
inst = odi.instrument(images[0])
#source = 'twomass'

#Create offlines catalogs
for img in images:
    print 'Retrieving QR SDSS and 2MASS catalogs for:', img
    for key in odi.OTA_dictionary:
        ota = odi.OTA_dictionary[key]
        outputsd = odi.sdsspath+'offline_'+ota+'.'+str(img[16:-5])+'.sdss'
        if not os.path.isfile(outputsd):
            x,y = odi.get_sdss_coords_offline(img,ota,inst,output=outputsd)
            output2m = odi.twomasspath+'offline_'+ota+'.'+str(img[16:-5])+'.mass'
            x,y = odi.get_2mass_coords_offline(img,ota,inst,output=output2m)

listfiles = glob.glob('*.lis')
if len(listfiles) == 0:
    odi.imcombine_lists(images, filters)
else:
    print 'imcombine lists done'

if not os.path.isfile('bpms.done'):
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


if not os.path.isfile('derived_props.txt'):
    f1 = open('derived_props.txt','w+')
    print >> f1, '# img  ota  filter fwhm  zp_med  zp_std  bg_mean  bg_med  bg_std'
    for img in images:
        for key in tqdm(odi.OTA_dictionary):
            ota = odi.OTA_dictionary[key]
            hdulist = odi.fits.open(img)
            hdr = hdulist[0].header
            filt = hdr['filter']
            image_to_correct = img+'['+ota+']'
            correction_image = ota+'.'+filt+'.med.fits'
            corrected_image = 'illcor_'+ota+'.'+str(img[16:])
            if not os.path.isfile(odi.illcorpath+corrected_image):
                odi.illumination_corrections(image_to_correct, correction_image, corrected_image)
            gaps = odi.get_gaps(img, ota)
            reprojed_image = 'reproj_'+ota+'.'+str(img[16:])
            if not os.path.isfile(odi.reprojpath+reprojed_image):
                pixcrd3 = odi.list_wcs_coords(img, ota, gaps, inst,output=img[:-5]+'.'+ota+'.radec.coo', gmaglim=23., stars_only=True, offline = True, source = source)
                try:
                    odi.fix_wcs(img, ota, coords=img[:-5]+'.'+ota+'.radec.coo', iters=3)
                except:
                    print 'msccmatch failed, wait a second and try again'
                    time.sleep(1.0)
                    odi.fix_wcs(img, ota, coords=img[:-5]+'.'+ota+'.radec.coo', iters=3)
                odi.reproject_ota(img, ota, rad, decd)
            gaps = odi.get_gaps_rep(img, ota)
            odi.refetch_sdss_coords(img, ota, gaps, inst,gmaglim=21.5,offline = True,source=source)
            #run an additional refetch to get the xy for 2mass so they can be used for scaling
            odi.repoxy_offline(img, ota, gaps, inst,gmaglim=21.5,source='twomass')
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
else:
    imgnum,fwhm,zp_med, zp_std, bg_mean, bg_median, bg_std = np.loadtxt('derived_props.txt',usecols=(0,3,4,5,6,7,8),unpack=True)
    ota_d, filt_d = np.loadtxt('derived_props.txt',usecols=(1,2),unpack=True,dtype=str)
    finished = zip(imgnum,ota_d,filt_d)
    f1 = open('derived_props.txt','a+')
    for img in images:
        for key in tqdm(odi.OTA_dictionary):
            ota = odi.OTA_dictionary[key]
            hdulist = odi.fits.open(img)
            hdr = hdulist[0].header
            filt = hdr['filter']
            finishcheck = (int(str(img[16])),ota,filt)
            if finishcheck in finished:
                already = 0
            else:
                image_to_correct = img+'['+ota+']'
                correction_image = ota+'.'+filt+'.med.fits'
                corrected_image = 'illcor_'+ota+'.'+str(img[16:])
                if not os.path.isfile(odi.illcorpath+corrected_image):
                    odi.illumination_corrections(image_to_correct, correction_image, corrected_image)
                gaps = odi.get_gaps(img, ota)
                reprojed_image = 'reproj_'+ota+'.'+str(img[16:])
                if not os.path.isfile(odi.reprojpath+reprojed_image):
                    pixcrd3 = odi.list_wcs_coords(img, ota, gaps, inst,output=img[:-5]+'.'+ota+'.radec.coo', gmaglim=23., stars_only=True, offline = True, source = source)
                    try:
                        odi.fix_wcs(img, ota, coords=img[:-5]+'.'+ota+'.radec.coo', iters=1)
                    except IndexError:     
                        print 'not enough stars to fix wcs, skipping for this ota:', img, ota
                    except:
                        print 'msccmatch failed, wait a second and try again'
                        time.sleep(1.0)
                        odi.fix_wcs(img, ota, coords=img[:-5]+'.'+ota+'.radec.coo', iters=1)
                    
                    odi.reproject_ota(img, ota, rad, decd)
                gaps = odi.get_gaps_rep(img, ota)
                odi.refetch_sdss_coords(img, ota, gaps, inst,gmaglim=21.5,offline = True,source=source)
                #run an additional refetch to get the xy for 2mass so they can be used for scaling
                odi.repoxy_offline(img, ota, gaps, inst,gmaglim=21.5,source='twomass')
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
