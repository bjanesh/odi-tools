import sys, os, glob, string
import numpy as np
import astropy as ast
import photutils as phot
import reproject as rep
import matplotlib.pyplot as plt
from astropy.visualization.mpl_normalize import ImageNormalize
import odi_calibrate as odi
from pyraf import iraf
from scipy.ndimage import binary_dilation
from astropy.io import fits
from astropy.wcs import WCS
#from rand_bkg import bkg_boxes
from astropy.convolution import Gaussian2DKernel
from astropy.stats import sigma_clipped_stats
from photutils.detection import detect_sources

OTA_dictionary = {1:'OTA33.SCI',2: 'OTA34.SCI',3 :'OTA44.SCI', 4:'OTA43.SCI',5:'OTA42.SCI', 6:'OTA32.SCI', 
          7:'OTA22.SCI' ,8:'OTA23.SCI',9:'OTA24.SCI'}

def ota_sdss(img,ota):
  full_hdulist = fits.open(img.f)
  sdss_cat = full_hdulist['CAT.PHOTCALIB']
  ota_cat = sdss_cat.data[np.where(sdss_cat.data['ODI_OTA'] == float(ota))]
  full_hdulist.close()
  return ota_cat

# Example of how to use if you want to create a new
# dictionary containing all otas.
img = '20140406T214040.1_GCPair-F1_odi_g.5914.fits'
cats = {}
for key in OTA_dictionary:
  ota_string = OTA_dictionary[key]
  ota_string = ota_string.strip('OTA.SCI')
  cats['sdss_'+ota_string] = ota_sdss(img,ota_string)

# The data then can be acceses with something like
# cats['sdss_33']['column_name']