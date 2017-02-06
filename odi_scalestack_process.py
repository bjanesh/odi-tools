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
    object_str, filters, instrument, images, illcor_flag, skyflat_src, wcs_flag, reproject_flag, scale_flag, stack_flag, gaia_flag, cluster_flag, ra_center, dec_center, min_radius = odi.cfgparse('config.yaml')
except IOError:
    print 'config.yaml does not exist, quitting...'
    exit()
    
source = 'sdss'
inst = odi.instrument(images[filters[0]][1])

imgnum,fwhm_d,zp_med, zp_std, bg_mean, bg_median, bg_std = np.loadtxt('derived_props.txt',usecols=(0,3,4,5,6,7,8),unpack=True)
ota_d, filt_d = np.loadtxt('derived_props.txt',usecols=(1,2),unpack=True,dtype=str)
id_d = zip(imgnum,ota_d,filt_d)
fwhm_dict = dict(zip(id_d,fwhm_d))

for filter in filters:
    # Scaling with all sources
    images_ = images[filter]
    print 'Scaling images for filter ',filter
    for img in images_:
        # img = images_[dith]
        dither  = img.dither()+'_'
        for key in tqdm(odi.OTA_dictionary):
            ota = odi.OTA_dictionary[key]
            if not os.path.isfile(odi.sourcepath+'source_'+ota+'.'+img.base()+'.csv'):
                odi.source_find(img,ota,inst)
                gaps = odi.get_gaps_rep(img, ota)
                odi.source_xy(img,ota,gaps,filter,inst)
                fwhm = odi.getfwhm_source(img,ota)
                #fwhm = fwhm_dict[img_id]
                print fwhm
                odi.phot_sources(img, ota, fwhm)
                odi.phot_combine(img, ota)
        if not os.path.isfile(odi.sourcepath+dither+filter+'.allsource'):
            dither_total = odi.sourcepath+dither+filter+'.allsource' 
            cat_command = 'cat' + ' ' + 'sources/'+'*'+dither+'*_'+filter+'*.totphot' + '>'+dither_total
            os.system(cat_command)
    
    # choose the initial reference image (lowest airmass to start)
    # print images_.values()
    refimg_ = odi.find_ref_image(images_)
    ref_img = images_[refimg_+1]
    
    # calculate scaling factors
    scales_ = {}
    stds_ = {}
    n_ = {}
    for dith in images_:
        img = images_[dith]
        scale,std,n = odi.source_scale(img,ref_img,filter)
        scales_[img] = scale
        stds_[img] = std
        n_[img] = n
    
    # recalculate scaling factors IF the highest scaling factor is not the initial reference image
    # print the scaling factors out to a file for review

    ims = scales_.keys()
    scls = scales_.values()
    new_ref = ims[np.argmax(scls)]
    if new_ref != ref_img:  
        ref_img = new_ref  
        for dith in images_:
            img = images_[dith]
            scale,std,n = odi.source_scale(img,ref_img,filter)
            scales_[img] = scale
            stds_[img] = std
            n_[img] = n
                
    with open(filter+'_scales.txt','w+') as sclfile:
        print >> sclfile, 'image'+' '*(len(images_[1])-4)+'scale   std     n'
        for dith in images_:
            img = images_[dith]
            print >> sclfile, img, '{0:7.5f} {1:7.5f} {2:5d}'.format(scales_[img], stds_[img], n_[img])
    
    # actually apply the scaling factors to the images
    if scale_flag:    
        for dith in images_:
            img = images_[dith]
            for key in odi.OTA_dictionary:
                ota = odi.OTA_dictionary[key]
                if not os.path.isfile(odi.scaledpath+'scaled_'+ota+'.'+img.stem()):
                    # gaps = odi.get_gaps_rep(img, ota)
                    odi.scale_ota(img, ota, scales_[img])
                    odi.force_update_bpm(img, ota)
    else:
        print 'scaling not performed, set flag in config.yaml'
    
    # finally stack the images
    if stack_flag:
        stacked_img = odi.stack_images(object_str, ref_img)
    else:
        print 'stacking not performed, set flag in config.yaml'
