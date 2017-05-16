import sys, os, glob, string
import numpy as np
import astropy as ast
import matplotlib.pyplot as plt
from pyraf import iraf
from tqdm import tqdm

import odi_config as odi

def dark_sky_flat(filter):
    # from cv2 import medianBlur
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
        
    iraf.set(clobber = 'yes')
    print 'smoothing dark sky flats for',filter
    for key in tqdm(odi.OTA_dictionary):
        image_list = odi.OTA_dictionary[key]+'.'+filter+'.lis'
        med_out = image_list.replace('.lis','.med.fits')
        # don't overwrite the original flat
        med_smooth = image_list.replace('.lis','.med.smooth.fits')

        # smooth the flat with a 50 x 50 median box
        iraf.unlearn(iraf.imutil.imarith,iraf.imfilter.median)
        iraf.imfilter.median.setParam('input',odi.skyflatpath+med_out)
        iraf.imfilter.median.setParam('output',odi.skyflatpath+med_smooth)
        iraf.imfilter.median.setParam('xwindow',51)
        iraf.imfilter.median.setParam('ywindow',51)
        iraf.imfilter.median.setParam('zloreject',1.0) # ignore 0.0s in the smoothing, they'll cause image artifacts
        iraf.imfilter.median(verbose='no', mode='h')
        
        # determine the normalization factor from the _smoothed_ image
        if key == 1:
            data,header = odi.fits.getdata(odi.skyflatpath+med_smooth,header=True)
            mean, median, std = odi.sigma_clipped_stats(data, sigma=3.0, mask_value=0.0) # be sure to ignore 0.0s in the flat
            normalization_factor = median
        
        # smoothing using numpy, this method is much slower
        # data, header = odi.fits.getdata(odi.skyflatpath+med_out, header=True)
        # # data_gaps = data.astype(bool)
        # data_smoothed = medianBlur(data, ksize=51)
        # plt.imshow(data_smoothed, interpolation='nearest')
        # plt.show()
        
        # divide the smoothed image by the normalization factor
        iraf.imutil.imarith.setParam('operand1',odi.skyflatpath+med_smooth)
        iraf.imutil.imarith.setParam('op','/')
        iraf.imutil.imarith.setParam('operand2',normalization_factor)
        iraf.imutil.imarith.setParam('result',odi.skyflatpath+med_smooth)
        iraf.imutil.imarith(verbose='no', mode='h')
        
        # make an image of 0s and 1s to restore the cell gaps
        # iraf.imutil.imexpr('(a == 0) ? 0 : 1', odi.skyflatpath+'temp_gaps', odi.skyflatpath+med_smooth, verbose='no')
        # iraf.imutil.imarith.setParam('operand1',odi.skyflatpath+med_smooth)
        # iraf.imutil.imarith.setParam('op','*')
        # iraf.imutil.imarith.setParam('operand2',odi.skyflatpath+'temp_gaps')
        # iraf.imutil.imarith.setParam('result',odi.skyflatpath+med_smooth)
        # iraf.imutil.imarith(verbose='no', mode='h')
        iraf.imutil.imdelete('temp_gaps', mode='h')
    iraf.set(clobber = 'no')

    return normalization_factor

def illumination_corrections(image_to_correct, correction_image, corrected_image,do_correction=True):
    #print image_to_correct,correction_image,corrected_image
    iraf.unlearn(iraf.imutil.imarith,iraf.imfilter.median)
    iraf.imutil.imarith.setParam('operand1',image_to_correct)
    iraf.imutil.imarith.setParam('op','/')
    if do_correction == True:
        iraf.imutil.imarith.setParam('operand2',odi.skyflatpath+correction_image)
    else:
        print 'not applying illcor'
        iraf.imutil.imarith.setParam('operand2',1.0)
    iraf.imutil.imarith.setParam('result',odi.illcorpath+corrected_image)
    iraf.imutil.imarith(mode='h')

def main():
    odi.OTA_dictionary = odi.podi_dictionary
    norm = dark_sky_flat('odi_g')
    print norm

if __name__ == '__main__':
    main()
