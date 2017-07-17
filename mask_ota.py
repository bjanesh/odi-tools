import sys, os, glob, string
import numpy as np
import astropy as ast
import matplotlib.pyplot as plt
from pyraf import iraf
import odi_config as odi


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
          hdu.writeto(odi.bppath+mask_name,clobber=True)

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
    pass

if __name__ == '__main__':
    main()
