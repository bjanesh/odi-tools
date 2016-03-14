from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.io import fits
from astropy.wcs import WCS
#from rand_bkg import bkg_boxes
from astropy.convolution import Gaussian2DKernel
from astropy.stats import sigma_clipped_stats
from photutils.detection import detect_sources
from photutils.utils import random_cmap
from scipy.ndimage import binary_dilation
import pandas as pd
import numpy as np
import os
import glob
import odi_config as odi

def get_sdss_coords_offline(img, ota, inst,output='test.sdss'):
    hdulist = odi.fits.open(img)
    hdu = hdulist[ota]
    
    if inst == 'podi':
        pvlist = hdu.header['PV*']
        for pv in pvlist:
            tpv = 'T'+pv
            hdu.header.rename_keyword(pv, tpv, force=False)
    xdim = hdu.header['NAXIS1']
    ydim = hdu.header['NAXIS2']
    
    sdss_cat_img = hdulist['CAT.PHOTCALIB']
    sdss_cat_img_df = pd.DataFrame.from_dict(sdss_cat_img.data)
    hdulist.close()
    
    ota = float(ota.strip('OTA.SCI'))
    ota_matches_df = sdss_cat_img_df.iloc[np.where(sdss_cat_img_df['ODI_OTA'] == ota)]

    needed_columns = ['SDSS_RA','SDSS_DEC','SDSS_MAG_U',
		      'SDSS_ERR_U', u'SDSS_MAG_G', u'SDSS_ERR_G', u'SDSS_MAG_R',
		      'SDSS_ERR_R', u'SDSS_MAG_I', u'SDSS_ERR_I', u'SDSS_MAG_Z',
		      'SDSS_ERR_Z','ODI_OTA']

    output_df = ota_matches_df[needed_columns]
    output_df.to_csv(output,index=False)
    return xdim, ydim

def get_2mass_coords_offline(img, ota, inst,output='test.mass'):
    hdulist = odi.fits.open(img)
    hdu = hdulist[ota]
    
    if inst == 'podi':
        pvlist = hdu.header['PV*']
        for pv in pvlist:
            tpv = 'T'+pv
            hdu.header.rename_keyword(pv, tpv, force=False)
    xdim = hdu.header['NAXIS1']
    ydim = hdu.header['NAXIS2']
    
    mass_cat_img = hdulist['CAT.ODI+2MASS']
    mass_cat_img_df = mass_cat_img.data
    hdulist.close()
        
    ota = float(ota.strip('OTA.SCI'))
    ota_matches_df = mass_cat_img_df[np.where(mass_cat_img_df['OTA'] == ota)]
    
    junk_df = pd.DataFrame.from_dict(ota_matches_df)
    
    output_df = junk_df
    output_df.to_csv(output,index=False)
    return xdim, ydim


def main():
    pass

if __name__ == '__main__':
    main()
