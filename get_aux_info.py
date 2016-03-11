#! /usr/bin/env python
from astropy.io import fits as pyfits
from astropy.table import Table, Column
import numpy as np
import subprocess
import os
import glob
import tarfile
from collections import OrderedDict


object_name = 'm13-se'
tar_list = glob.glob('*_'+object_name+'_*.tar')
tar_list.sort()



frame             = []
airmass          =  []
exptime          =  []
filter           =  []
fwhm_sdss_arcsec =  []
fwhm_sdss_pix    =  []
fwhm_all_arcsec  =  []
fwhm_all_pix     =  []
skynoise         =  []
skylevel         =  []                                                 
skybgstd         =  []


for i in range(len(tar_list)):
  tar = tarfile.open(tar_list[i])
  tar.extractall()
  tar.close()
  aux_path = tar_list[i].strip('.tar')
  print aux_path
  os.chdir(aux_path)
  os.chdir('auxiliary_products')
  current_frame = glob.glob('*.fits.fz')
  hdulist = pyfits.open(current_frame[0])
  header           =  hdulist[0].header
  frame.append(current_frame[0])
  airmass.append(header['airmass'])
  exptime.append(header['exptime'])
  filter.append(header['filter'])
  fwhm_sdss_arcsec.append(header['fwhmstar'])
  fwhm_sdss_pix.append(fwhm_sdss_arcsec[i]/0.11)
  fwhm_all_arcsec.append(header['fwhm_all'])
  fwhm_all_pix.append(fwhm_all_arcsec[i]/0.11)
  skynoise.append(header['skynoise'])
  skylevel.append(header['skylevel'])                                                    
  skybgstd.append(header['skybgstd'])
  os.chdir('../..')
  
aux_data = OrderedDict([('frame',frame),('filter',filter),('exptime',exptime),('airmass',airmass),
			('fwhm_all_arcsec',fwhm_all_arcsec),('fwhm_sdss_arcsec',fwhm_sdss_arcsec),('fwhm_all_pix',fwhm_all_pix),
			('fwhm_sdss_pix',fwhm_sdss_pix),('skynoise',skynoise),('skylevel',skylevel),
			('skybgstd',skybgstd)])
  

  
aux_table = Table(aux_data)
aux_table.write(object_name+'.aux',format='ascii.commented_header')