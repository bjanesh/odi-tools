import sys, os, glob, string
import numpy as np
import astropy as ast
import matplotlib.pyplot as plt
from pyraf import iraf
from tqdm import tqdm

import odi_config as odi

def dark_sky_flat(filter):
    med_otalist = []
    print 'making dark sky flats for',filter
    for key in tqdm(odi.OTA_dictionary):
        image_list = odi.OTA_dictionary[key]+'.'+filter+'.lis'
        med_out = image_list.replace('.lis','.med.fits')
        med_otalist.append(med_out)
        iraf.unlearn(iraf.immatch.imcombine)
        iraf.immatch.imcombine.setParam('input','@'+str(image_list))
        iraf.immatch.imcombine.setParam('output',odi.skyflatpath+med_out)
        iraf.immatch.imcombine.setParam('combine','median')
        iraf.immatch.imcombine.setParam('masktype','goodvalue')
        iraf.immatch.imcombine.setParam('maskvalue',0)
        iraf.immatch.imcombine.setParam('scale','median')
        # iraf.immatch.imcombine.setParam('zero','none')
        iraf.immatch.imcombine.setParam('zero','median')
        iraf.immatch.imcombine(logfile='imcombine.log.txt', mode='h')
        if key == 1:
            data,header = odi.fits.getdata(odi.skyflatpath+med_out,header=True)
            mean, median, std = odi.sigma_clipped_stats(data, sigma=3.0)
            normalization_factor = median
    iraf.set(clobber = 'yes')
    print 'smoothing dark sky flats for',filter
    for i in tqdm(range(len(med_otalist))):
        iraf.unlearn(iraf.imutil.imarith,iraf.imfilter.median)
        iraf.imutil.imarith.setParam('operand1',odi.skyflatpath+med_otalist[i])
        iraf.imutil.imarith.setParam('op','/')
        iraf.imutil.imarith.setParam('operand2',normalization_factor)
        iraf.imutil.imarith.setParam('result',odi.skyflatpath+med_otalist[i])
        iraf.imutil.imarith(verbose='no', mode='h')
        iraf.imfilter.median.setParam('input',odi.skyflatpath+med_otalist[i])
        iraf.imfilter.median.setParam('output',odi.skyflatpath+med_otalist[i])
        iraf.imfilter.median.setParam('xwindow',3)
        iraf.imfilter.median.setParam('ywindow',3)
        iraf.imfilter.median(verbose='no', mode='h')
    iraf.set(clobber = 'no')

    return normalization_factor

def illumination_corrections(image_to_correct, correction_image, corrected_image):
    #print image_to_correct,correction_image,corrected_image
    iraf.unlearn(iraf.imutil.imarith,iraf.imfilter.median)
    iraf.imutil.imarith.setParam('operand1',image_to_correct)
    iraf.imutil.imarith.setParam('op','/')
    iraf.imutil.imarith.setParam('operand2',odi.skyflatpath+correction_image)
    iraf.imutil.imarith.setParam('result',odi.illcorpath+corrected_image)
    iraf.imutil.imarith(mode='h')

def main():
	pass

if __name__ == '__main__':
	main()

