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

try:
    object_str, filters, instrument, images, illcor_flag, skyflat_src, wcs_flag, reproject_flag, scale_flag, stack_flag, gaia_flag, cluster_flag, ra_center, dec_center, min_radius = odi.cfgparse('config.yaml')
except IOError:
    print 'config.yaml does not exist, quitting...'
    exit()

# for basic processing, filter shouldn't matter enough to consider separately
# or rather, this script already handles that just fine
# so just stick all the image names together into one long list
images_ = [img for sublist in images.values() for img in sublist]
# print images_

rad, decd = odi.get_targ_ra_dec(images_[0], 'OTA33.SCI')
if gaia_flag:
    source = 'gaia'
else:
    source = 'sdss'
inst = odi.instrument(images_[0])
#source = 'twomass'

#Create offline catalogs
for img in images_:
    print 'Retrieving QR SDSS and Gaia catalogs for:', img.stem()
    for key in tqdm(odi.OTA_dictionary):
        ota = odi.OTA_dictionary[key]
        outputsd = odi.sdsspath+'offline_'+ota+'.'+img.base()+'.sdss'
        if not os.path.isfile(outputsd):
            x,y = odi.get_sdss_coords_offline(img,ota,inst,output=outputsd)
            output2m = odi.twomasspath+'offline_'+ota+'.'+img.base()+'.mass'
            # x,y = odi.get_2mass_coords_offline(img,ota,inst,output=output2m)
        if gaia_flag == True:
            outputg = odi.gaiapath+'offline_'+ota+'.'+img.base()+'.gaia'
            if not os.path.isfile(outputg):
                if cluster_flag == True:
                    odi.get_gaia_coords(img,ota,inst,
                                        output=outputg,
                                        cluster=cluster_flag,
                                        racenter=float(ra_center),
                                        deccenter=float(dec_center),
                                        min_radius=float(min_radius),
                                        G_lim = 20.75)
                else:
                    odi.get_gaia_coords(img,ota,inst,
                                        output=outputg,
                                        cluster=cluster_flag)



listfiles = glob.glob('*.lis')
if len(listfiles) == 0:
    odi.imcombine_lists(images_, filters)
else:
    print 'imcombine lists done'

if not os.path.isfile('bpms.done'):
    for img in images_:
        print 'updating bpms for', img.stem()
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
    for img in images_:
        for key in tqdm(odi.OTA_dictionary):
            ota = odi.OTA_dictionary[key]
            hdulist = odi.fits.open(img.f)
            hdr = hdulist[0].header
            filt = hdr['filter']
            image_to_correct = img.f+'['+ota+']'
            correction_image = ota+'.'+filt+'.med.fits'
            corrected_image = 'illcor_'+ota+'.'+img.stem()
            if not os.path.isfile(odi.illcorpath+corrected_image):
                odi.illumination_corrections(image_to_correct, correction_image, corrected_image, do_correction=illcor_flag)
            gaps = odi.get_gaps(img, ota)
            reprojed_image = 'reproj_'+ota+'.'+img.stem()
            if not os.path.isfile(odi.reprojpath+reprojed_image):
                if wcs_flag:
                    pixcrd3 = odi.list_wcs_coords(img, ota, gaps, inst,output=img.nofits()+'.'+ota+'.radec.coo', gmaglim=23., stars_only=True, offline = True, source = source)
                    try:
                        odi.fix_wcs(img, ota, coords=img.nofits()+'.'+ota+'.radec.coo', iters=1)
                    except:
                        try:
                            print 'msccmatch failed, wait a second and try again'
                            time.sleep(1.0)
                            odi.fix_wcs(img, ota, coords=img.nofits()+'.'+ota+'.radec.coo', iters=1)
                        except:
                            print 'there might be too few stars for msccmatch, just skip it.'
                if reproject_flag:
                    wcsref = odi.illcorpath+'illcor_OTA33.SCI.'+images_[0].stem()
                    odi.reproject_ota(img, ota, rad, decd, wcsref)
            gaps = odi.get_gaps_rep(img, ota)
            odi.refetch_sdss_coords(img, ota, gaps, inst,gmaglim=21.5,offline = True,source=source)
            #run an additional refetch to get the xy for 2mass so they can be used for scaling
            # odi.repoxy_offline(img, ota, gaps, inst,gmaglim=21.5,source='twomass')
            odi.repoxy_offline(img, ota, gaps, inst,gmaglim=21.5,source='gaia')
            fwhm = odi.getfwhm_ota(img, ota, gaia=gaia_flag)
            if 'odi_NB695' in filters:
                zp_med, zp_std = 99.99,99.99
            elif source == 'sdss':
                zp_med, zp_std = 99.99,99.99
            elif source == 'twomass':
                zp_med, zp_std = 99.99,99.99
            elif source == 'gaia':
                zp_med, zp_std = 99.99,99.99
            if not os.path.isfile(odi.bgsubpath+'bgsub_'+ota+'.'+img.stem()):
                bg_mean, bg_median, bg_std = odi.bgsub_ota(img, ota, apply=True)
            else:
                bg_mean, bg_median, bg_std = odi.bgsub_ota(img, ota, apply=False)
            print >> f1, img.dither(), ota, filt, fwhm, zp_med, zp_std, bg_mean, bg_median, bg_std
            dim_stats = odi.check_mask_dim(img,ota)
            if not dim_stats:
                print 'mask dimensions do not match image'
                print 'redo', img, ota
                raise ValueError
    f1.close()
else:
    imgnum,fwhm,zp_med, zp_std, bg_mean, bg_median, bg_std = np.loadtxt('derived_props.txt',usecols=(0,3,4,5,6,7,8),unpack=True)
    ota_d, filt_d = np.loadtxt('derived_props.txt',usecols=(1,2),unpack=True,dtype=str)
    finished = zip(imgnum,ota_d,filt_d)
    f1 = open('derived_props.txt','a+')
    for img in images_:
        for key in tqdm(odi.OTA_dictionary):
            ota = odi.OTA_dictionary[key]
            hdulist = odi.fits.open(img.f)
            hdr = hdulist[0].header
            filt = hdr['filter']
            finishcheck = (int(img.dither()),ota,filt)
            if finishcheck in finished:
                already = 0
            else:
                image_to_correct = img.f+'['+ota+']'
                correction_image = ota+'.'+filt+'.med.fits'
                corrected_image = 'illcor_'+ota+'.'+img.stem()
                if not os.path.isfile(odi.illcorpath+corrected_image):
                    odi.illumination_corrections(image_to_correct, correction_image, corrected_image, do_correction=illcor_flag)
                gaps = odi.get_gaps(img, ota)
                reprojed_image = 'reproj_'+ota+'.'+img.stem()
                if not os.path.isfile(odi.reprojpath+reprojed_image):
                    if wcs_flag:
                        pixcrd3 = odi.list_wcs_coords(img, ota, gaps, inst,output=img.nofits()+'.'+ota+'.radec.coo', gmaglim=23., stars_only=True, offline = True, source = source)
                        try:
                            odi.fix_wcs(img, ota, coords=img.nofits()+'.'+ota+'.radec.coo', iters=1)
                        except IndexError:
                            print 'not enough stars to fix wcs, skipping for this ota:', img, ota
                        except:
                            print 'msccmatch failed, wait a second and try again'
                            time.sleep(1.0)
                            odi.fix_wcs(img, ota, coords=img.nofits()+'.'+ota+'.radec.coo', iters=1)
                    if reproject_flag:
                        wcsref = odi.illcorpath+'illcor_OTA33.SCI.'+images_[0].stem()
                        odi.reproject_ota(img, ota, rad, decd, wcsref)
                gaps = odi.get_gaps_rep(img, ota)
                odi.refetch_sdss_coords(img, ota, gaps, inst,gmaglim=21.5,offline = True,source=source)
                #run an additional refetch to get the xy for 2mass so they can be used for scaling
                # odi.repoxy_offline(img, ota, gaps, inst,gmaglim=21.5,source='twomass')
                odi.repoxy_offline(img, ota, gaps, inst,gmaglim=21.5,source='gaia')
                fwhm = odi.getfwhm_ota(img, ota, gaia=gaia_flag)
                if source == 'sdss':
                    zp_med, zp_std = 99.99,99.99
                if source == 'twomass':
                    zp_med, zp_std = 99.99,99.99
                elif source == 'gaia':
                    zp_med, zp_std = 99.99,99.99
                if not os.path.isfile(odi.bgsubpath+'bgsub_'+ota+'.'+img.stem()):
                    bg_mean, bg_median, bg_std = odi.bgsub_ota(img, ota, apply=True)
                else:
                    bg_mean, bg_median, bg_std = odi.bgsub_ota(img, ota, apply=False)
                print >> f1, img.dither(), ota, filt, fwhm, zp_med, zp_std, bg_mean, bg_median, bg_std
                dim_stats = odi.check_mask_dim(img,ota)
                if not dim_stats:
                    print 'mask dimensions do not match image'
                    print 'redo', img.stem(), ota
                    raise ValueError
    f1.close()
