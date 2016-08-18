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

images = glob.glob('*.fits')
filters = ['odi_g','odi_r','odi_i']

inst = odi.instrument(images[0])
#source = 'twomass'

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

for filt in filters:
    med_otalist = []
    print 'making master dark sky flats for',filt
    for key in tqdm(odi.OTA_dictionary):
        image_list = odi.OTA_dictionary[key]+'.'+filt+'.lis'
        med_out = image_list.replace('.lis','.mastersflat.fits')
        med_otalist.append(med_out)
        iraf.unlearn(iraf.immatch.imcombine)
        iraf.immatch.imcombine.setParam('input','@'+str(image_list))
        iraf.immatch.imcombine.setParam('output',med_out)
        iraf.immatch.imcombine.setParam('combine','median')
        iraf.immatch.imcombine.setParam('masktype','goodvalue')
        iraf.immatch.imcombine.setParam('maskvalue',0)
        iraf.immatch.imcombine.setParam('scale','median')
        # iraf.immatch.imcombine.setParam('zero','none')
        iraf.immatch.imcombine.setParam('zero','median')
        iraf.immatch.imcombine(logfile='imcombine.log.txt', mode='h')
        if key == 1:
            data,header = odi.fits.getdata(med_out,header=True)
            mean, median, std = odi.sigma_clipped_stats(data, sigma=3.0)
            normalization_factor = median
    iraf.set(clobber = 'yes')
    print 'smoothing dark sky flats for',filt
    for i in tqdm(range(len(med_otalist))):
        iraf.unlearn(iraf.imutil.imarith,iraf.imfilter.median)
        iraf.imutil.imarith.setParam('operand1',med_otalist[i])
        iraf.imutil.imarith.setParam('op','/')
        iraf.imutil.imarith.setParam('operand2',normalization_factor)
        iraf.imutil.imarith.setParam('result',med_otalist[i])
        iraf.imutil.imarith(verbose='no', mode='h')
        iraf.imfilter.median.setParam('input',med_otalist[i])
        iraf.imfilter.median.setParam('output',med_otalist[i])
        iraf.imfilter.median.setParam('xwindow',3)
        iraf.imfilter.median.setParam('ywindow',3)
        iraf.imfilter.median(verbose='no', mode='h')
iraf.set(clobber = 'no')