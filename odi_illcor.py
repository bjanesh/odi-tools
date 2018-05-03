import sys, os, glob, string
import numpy as np
from pyraf import iraf
from tqdm import tqdm

import odi_config as odi

def dark_sky_flat(filter, box_size=51):
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
        iraf.immatch.imcombine.setParam('zero','none')
        # iraf.immatch.imcombine.setParam('zero','median')
        iraf.immatch.imcombine(logfile='imcombine.log.txt', mode='h')
        
    iraf.set(clobber = 'yes')
    print 'smoothing dark sky flats for',filter,'with box size {0:3d}x{0:3d}'.format(box_size)
    for key in tqdm(odi.OTA_dictionary):
        image_list = odi.OTA_dictionary[key]+'.'+filter+'.lis'
        med_out = image_list.replace('.lis','.med.fits')
        # don't overwrite the original flat
        med_smooth = image_list.replace('.lis','.med.smooth.fits')

        # smooth the flat with a N x N median box
        iraf.unlearn(iraf.imutil.imarith,iraf.imfilter.median)
        iraf.imfilter.fmedian.setParam('input',odi.skyflatpath+med_out)
        iraf.imfilter.fmedian.setParam('output',repr(key)+'temp_smooth.fits')
        iraf.imfilter.fmedian.setParam('xwindow',box_size)
        iraf.imfilter.fmedian.setParam('ywindow',box_size)
        iraf.imfilter.fmedian.setParam('zloreject',1.0) # ignore 0.0s in the smoothing, they'll cause image artifacts
        iraf.imfilter.fmedian(verbose='no', mode='h')
        
        # smoothing is not leaving zeros but instead some very small number--replace them with 0.0s
        iraf.imutil.imexpr('(a < 1.51) ? 0 : a',odi.skyflatpath+med_smooth,repr(key)+'temp_smooth.fits',verbose='no')
        iraf.imutil.imdelete(repr(key)+'temp_smooth.fits')
        
        # determine the normalization factor from the _smoothed_ image
        if key == 1:
            data,header = odi.fits.getdata(odi.skyflatpath+med_smooth,header=True)
            mean, median, std = odi.sigma_clipped_stats(data, sigma=3.0, mask_value=0.0) # be sure to ignore 0.0s in the flat
            normalization_factor = median
            # print normalization_factor
        
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
        tqdm.write('not applying illcor')
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
            hdu.writeto(mask_name,overwrite=True)
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
    
def get_gaps(img, ota):
    """
    Create a numpy array mask of the gaps in an ota.

    Parameters
    ----------
    img : str
        Name of image
    ota : str
        Name of OTA
    Returns
    -------
    gaps_mask : numpy array
        A numpy array of the gap location on the ota.
    """
    hdulist = odi.fits.open(img.f)
    hdu = hdulist[ota]
    gaps_mask = (np.isnan(hdu.data)).astype(int)
    hdulist.close()
    return gaps_mask

def get_gaps_rep(img, ota):
    """
    Create a numpy array mask of the gaps in a reprojected ota.

    Parameters
    ----------
    img : str
        Name of image
    ota : str
        Name of OTA
    Returns
    -------
    gaps_mask : numpy array
        A numpy array of the gap location on the ota.
    """

    image = odi.reprojpath+'reproj_'+ota+'.'+img.stem()
    hdulist = odi.fits.open(image)
    hdu = hdulist[0]
    # plt.imshow(hdu.data, origin='lower', cmap='Greys_r', vmin=-10., vmax=500.)
    # plt.show()
    gaps_mask1 = (hdu.data<1.0).astype(int)
    selem = np.ones((5,5))    # dilate using a 25x25 box
    gaps_mask = odi.binary_dilation(gaps_mask1, selem)
    hdulist.close()

    # also update the bad pixel mask for the image to make sure the cell gaps are masked
    # this is necessary for the final imcombine
    mask_name = odi.bppath+'reproj_mask_'+ota+'.'+img.stem()
    BPM = mask_name.replace('fits','pl')
    if not os.path.isfile(BPM):
        # mask,gaps = mask_ota(img,ota)
        hdu = odi.fits.PrimaryHDU(gaps_mask.astype(float))
        if not os.path.isfile(mask_name):
            hdu.writeto(mask_name,overwrite=True)
    if not os.path.isfile(mask_name.replace('fits','pl')):
        iraf.unlearn(iraf.imutil.imcopy)
        iraf.imutil.imcopy.setParam('input',mask_name)
        iraf.imutil.imcopy.setParam('output',mask_name.replace('fits','pl'))
        iraf.imutil.imcopy.setParam('verbose','no')
        iraf.imutil.imcopy(mode='h')
    # if os.path.isfile(mask_name): # we don't need to keep the reproj fits mask, it takes up a ton of space
    #     iraf.imutil.imdelete(mask_name, verify='no', mode='h')
    iraf.unlearn(iraf.imutil.hedit)
    iraf.imutil.hedit.setParam('images',image)
    iraf.imutil.hedit.setParam('fields','BPM')
    iraf.imutil.hedit.setParam('value',BPM)
    iraf.imutil.hedit.setParam('add','yes')
    iraf.imutil.hedit.setParam('addonly','no')
    iraf.imutil.hedit.setParam('verify','no')
    iraf.imutil.hedit.setParam('update','yes')
    iraf.imutil.hedit(mode='h')

    return gaps_mask
    
def mask_ota(img, ota, reproj=False, deep_obj=False):
    """
    Create a numpy array that will be used to create the bad pixel mask
    for an ota. The function finds the pixel locations of the gaps in the
    ota as well as hot and dead pixels.

    Parameters
    -----------
    img : str
       String containing name of the image currently in use

    ota : str
       Name of ota extension to be used (e.g. OTA33.SCI)

    reproj : boolean
        If ``reproj`` is ``True`` this function returns background statistics on
        the OTA, but the mask is still produced.

    deep_obj: boolean
        If ``deep_obj`` is ``True`` the threshold for dead pixels is set to -900.

    Returns
    --------
    total_mask : 2D array
               Array mask for hot pixels, dead pixels, and the gaps.
    gap_mask : 2D array
             Array mask for the gaps.

    Notes
    -----
    Here are the limits used to located the pixels that are to be masked

    ``hot_pixels`` : hdu_ota.data > 58000.0

    ``dead_pixels`` : hdu_ota.data < 1.0

    ``gap_pixels`` : hdu_ota.data = ``NaN``


    """

    if reproj:
      image = odi.reprojpath+'reproj_'+ota+'.'+img.stem()
      QR_raw = odi.fits.open(image)
      hdu_ota = QR_raw[0]
    elif deep_obj:
      image = odi.otastackpath+ota+'_stack.fits'
      refimg = odi.scaledpath+'scaled_'+ota+'.'+img.stem()
      QR_raw = odi.fits.open(image)
      hdu_ota = QR_raw[0]
    else:
      QR_raw = odi.fits.open(img.f)
      hdu_ota = QR_raw[ota]

    # Mask hot pixels count greater than 58000
    hot_pix_mask = (hdu_ota.data > 58000.0).astype(int)
    # Mask dead pixels
    if deep_obj:
      dead_pix_mask = (hdu_ota.data < -900.0).astype(int)
    else:
      dead_pix_mask = (hdu_ota.data < 1.0).astype(int)
    # Mask gaps
    gaps_mask = np.isnan(hdu_ota.data).astype(int)

    #Calculate background stats
    bg_stats,bg_median,med_std,std_std,centers,max_box = odi.bkg_boxes(hdu_ota,100,20.0,sources=True)

    #print bg_median,med_std,std_std

    threshold = bg_median + (med_std * 3.0)
    segm_img = odi.detect_sources(hdu_ota.data, threshold, npixels=25)
    source_mask1 =  segm_img.data.astype(np.bool)
    selem = np.ones((15,15))    # dilate using a 25x25 box
    source_mask2 = odi.binary_dilation(source_mask1, selem)

    full_mask = source_mask2 + hot_pix_mask + dead_pix_mask + gaps_mask
    #full_mask = segm_img.data + hot_pix_mask + dead_pix_mask + gaps_mask
    total_mask = full_mask.astype(np.bool) # turn segm_img into a mask

    #create numpy masked array object
    # masked_array = ma.masked_array(hdu_ota.data,mask=total_mask)

    if reproj:
      # if operating on the reprojected image, return background statistics instead of the mask
      # but use the mask to get rid of sources, so it's a clean measurement
      mean, median, std = odi.sigma_clipped_stats(hdu_ota.data, mask=total_mask, sigma=3.0, iters=1)
      return mean, median, std#, dead_pix_mask
    if deep_obj:
      mask_name = 'objmask_'+ota+'.fits'
      # BPM = mask_name.replace('fits','pl')
      if not os.path.isfile(odi.bppath+mask_name):
          hdu = odi.fits.PrimaryHDU(source_mask2.astype(float))
          hdu.writeto(odi.bppath+mask_name,overwrite=True)

      ota_mask = 'objmask_'+ota+'.'+str(img.dither())+'.fits'

      if not os.path.isfile(odi.bppath+ota_mask):
          iraf.unlearn(iraf.mscred.mscimage)
          iraf.mscred.mscimage.format='image'
          iraf.mscred.mscimage.pixmask='no'
          iraf.mscred.mscimage.verbose='no'
          iraf.mscred.mscimage.wcssource='match'
          iraf.mscred.mscimage.reference=refimg
          # iraf.mscred.mscimage.ra=rad
          # iraf.mscred.mscimage.dec=decd
          iraf.mscred.mscimage.scale=0.11
          iraf.mscred.mscimage.rotation=0.0
          iraf.mscred.mscimage.blank=0
          iraf.mscred.mscimage.interpo='linear'
          # iraf.mscred.mscimage.minterp='poly5'
          iraf.mscred.mscimage.nxbl=2048
          iraf.mscred.mscimage.nybl=2048
          iraf.mscred.mscimage.fluxcon='yes'
          iraf.mscred.mscimage(odi.bppath+mask_name,odi.bppath+ota_mask)

      return
    else:
      return total_mask, gaps_mask
    QR_raw.close()
    

def main():
    odi.OTA_dictionary = odi.podi_dictionary
    norm = dark_sky_flat('odi_g')
    print norm

if __name__ == '__main__':
    main()
