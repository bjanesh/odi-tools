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
    
def force_update_bpm(img, ota):
    image = odi.scaledpath+'scaled_'+ota+'.'+img.stem()
    mask_name = odi.bppath+'reproj_mask_'+ota+'.'+img.stem()
    BPM = mask_name.replace('fits','pl')
    iraf.unlearn(iraf.imutil.hedit)
    iraf.imutil.hedit.setParam('images',image)
    iraf.imutil.hedit.setParam('fields','BPM')
    iraf.imutil.hedit.setParam('value',BPM)
    iraf.imutil.hedit.setParam('add','yes')
    iraf.imutil.hedit.setParam('addonly','no')
    iraf.imutil.hedit.setParam('verify','no')
    iraf.imutil.hedit.setParam('update','yes')
    iraf.imutil.hedit(show='no', mode='h')
    if os.path.isfile(mask_name):
        os.remove(mask_name)
    return

def make_bpms(img, ota):
    # for i in range(len(images)):
    #   for key in OTA_dictionary:
    # create string for mask fits name
    mask_name = odi.bppath+'mask_'+ota+'.'+img.stem()
    BPM = mask_name.replace('fits','pl')
    if not os.path.isfile(BPM):
        mask,gaps = odi.mask_ota(img,ota)
        hdu = odi.fits.PrimaryHDU(mask.astype(float))
        if not os.path.isfile(mask_name):
            hdu.writeto(mask_name,clobber=True)
        #if not os.path.isfile(mask_name.replace('fits','pl')):
            iraf.unlearn(iraf.imutil.imcopy)
            iraf.imutil.imcopy.setParam('input',mask_name)
            iraf.imutil.imcopy.setParam('output',mask_name.replace('fits','pl'))
            iraf.imutil.imcopy.setParam('verbose','no')
            iraf.imutil.imcopy(mode='h')
    iraf.unlearn(iraf.imutil.hedit)
    iraf.imutil.hedit.setParam('images',img.f+'['+ota+']')
    iraf.imutil.hedit.setParam('fields','BPM')
    iraf.imutil.hedit.setParam('value',BPM)
    iraf.imutil.hedit.setParam('add','yes')
    iraf.imutil.hedit.setParam('addonly','no')
    iraf.imutil.hedit.setParam('verify','no')
    iraf.imutil.hedit.setParam('update','yes')
    iraf.imutil.hedit(show='no', mode='h')
    if os.path.isfile(mask_name):
        os.remove(mask_name)
    return

def get_dims(header_string):
    dims = header_string[header_string.index('['):
                         header_string.index(']')+1].split(',')
    dims1 = int(dims[0].strip('['))
    dims2 = int(dims[1].strip(']'))
    return dims1,dims2


def check_mask_dim(img,ota):
    mask_name = odi.bppath+'reproj_mask_'+ota+'.'+img.stem()
    reproj_img = odi.reprojpath+'reproj_'+ota+'.'+img.stem()
    bgsub_img = odi.bgsubpath+'bgsub_'+ota+'.'+img.stem()

    mask_header = iraf.imhead(mask_name,Stdout=1)[0]
    mask_dim1, mask_dim2 = get_dims(mask_header)

    reproj_header = iraf.imhead(reproj_img,Stdout=1)[0]
    reproj_dim1, reproj_dim2 = get_dims(reproj_header)

    bgsub_header = iraf.imhead(bgsub_img,Stdout=1)[0]
    bgsub_dim1, bgsub_dim2 = get_dims(bgsub_header)

    dim_status = ((mask_dim1 == reproj_dim1 == bgsub_dim1) &
                  (mask_dim2 == reproj_dim2 == bgsub_dim2))
    if dim_status == False:
        tqdm.write('mask and image dimensions to not match for: ', img, ota)
        tqdm.write('mask: {:5d} {:5d}'.format(mask_dim1,mask_dim2))
        tqdm.write('reproj: {:5d} {:5d}'.format(reproj_dim1,reproj_dim2))
        tqdm.write('bgsub: {:5d} {:5d}'.format(bgsub_dim1,bgsub_dim2))
    else:
        tqdm.write('dimension test passed')
        tqdm.write('mask: {:5d} {:5d}'.format(mask_dim1,mask_dim2))
        tqdm.write('reproj: {:5d} {:5d}'.format(reproj_dim1,reproj_dim2))
        tqdm.write('bgsub: {:5d} {:5d}'.format(bgsub_dim1,bgsub_dim2))
    return dim_status

def main():
    odi.OTA_dictionary = odi.podi_dictionary
    norm = dark_sky_flat('odi_g')
    print norm

if __name__ == '__main__':
    main()
