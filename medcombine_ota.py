from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.visualization import *
from astropy.visualization.mpl_normalize import ImageNormalize
from rand_bkg import bkg_boxes
from astropy.convolution import Gaussian2DKernel
from astropy.stats import sigma_clipped_stats
from photutils.detection import detect_sources
from scipy.ndimage import binary_dilation
import numpy.ma as ma
from mask_ota import mask_ota
from pyraf import iraf
import os
import sys

images = ['20130510T002928.1_m13-9_odi_r.5869.fits', '20130510T002928.3_m13-9_odi_r.5869.fits', '20130510T002928.2_m13-9_odi_r.5869.fits']

OTA_dictionary = {1:'OTA33.SCI',2: 'OTA34.SCI',3 :'OTA44.SCI', 4:'OTA43.SCI',5:'OTA42.SCI', 6:'OTA32.SCI', 
		  7:'OTA22.SCI' ,8:'OTA23.SCI',9:'OTA24.SCI'}

bpmdirectory = 'bpmasks'
if not os.path.exists(bpmdirectory):
  print 'Creating director for bad pixel masks...'
  os.makedirs(bpmdirectory)
  
bppath = bpmdirectory+'/'

illcordirectory = 'illcor'
if not os.path.exists(illcordirectory):
  print 'Creating director for illumination corrected ota images...'
  os.makedirs(illcordirectory)
  
illcorpath = illcordirectory+'/'
  
def make_bpms(images):
  for i in range(len(images)):
    for key in OTA_dictionary:
      mask,masked_array = mask_ota(images[i],OTA_dictionary[key])
      hdu = fits.PrimaryHDU(mask.astype(float))
      # create string for mask fits name
      mask_name = bppath+'mask_'+OTA_dictionary[key]+'.'+str(images[i][-27:49])
      BPM = mask_name.replace('fits','pl')
      if not os.path.isfile(mask_name):
        hdu.writeto(mask_name,clobber=True)
      if not os.path.isfile(mask_name.replace('fits','pl')):
        iraf.unlearn(iraf.imutil.imcopy)
        iraf.imutil.imcopy.setParam('input',mask_name)
        iraf.imutil.imcopy.setParam('output',mask_name.replace('fits','pl'))
        iraf.imutil.imcopy.setParam('verbose','no')
        iraf.imutil.imcopy(mode='h')
      iraf.unlearn(iraf.imutil.hedit)
      iraf.imutil.hedit.setParam('images',images[i]+'['+str(key)+']')
      iraf.imutil.hedit.setParam('fields','BPM')
      iraf.imutil.hedit.setParam('value',BPM)
      iraf.imutil.hedit.setParam('add','no')
      iraf.imutil.hedit.setParam('addonly','no')
      iraf.imutil.hedit.setParam('verify','no')
      iraf.imutil.hedit.setParam('update','yes')
      iraf.imutil.hedit(mode='h')
  return

def imcombine_lists(images):
  for key in OTA_dictionary:
    list_name =  open(OTA_dictionary[key]+'.lis',"w")
    for i in range(len(images)):
      print >>list_name,images[i]+'['+str(key)+']'
    list_name.close()
  return
imcombine_lists(images)

def dark_sky_flat():
  med_otalist = []
  for key in OTA_dictionary:
    image_list = OTA_dictionary[key]+'.lis'
    med_out = image_list.replace('.lis','.med.fits')
    med_otalist.append(med_out)
    iraf.unlearn(iraf.immatch.imcombine)
    iraf.immatch.imcombine.setParam('input','@'+str(image_list))
    iraf.immatch.imcombine.setParam('output',med_out)
    iraf.immatch.imcombine.setParam('combine','median')
    iraf.immatch.imcombine.setParam('masktype','goodvalue')
    iraf.immatch.imcombine.setParam('maskvalue',0)
    iraf.immatch.imcombine.setParam('scale','median')
    iraf.immatch.imcombine.setParam('zero','median')
    iraf.immatch.imcombine(mode='h')
    if key == 1:
      data,header = fits.getdata(med_out,header=True)
      mean, median, std = sigma_clipped_stats(data, sigma=3.0)
      normalization_factor = median
  iraf.set(clobber = 'yes')
  for i in range(len(med_otalist)):
    iraf.unlearn(iraf.imutil.imarith,iraf.imfilter.median)
    iraf.imutil.imarith.setParam('operand1',med_otalist[i])
    iraf.imutil.imarith.setParam('op','/')
    iraf.imutil.imarith.setParam('operand2',normalization_factor)
    iraf.imutil.imarith.setParam('result',med_otalist[i])
    iraf.imutil.imarith(mode='h')
    iraf.imfilter.median.setParam('input',med_otalist[i])
    iraf.imfilter.median.setParam('output',med_otalist[i])
    iraf.imfilter.median.setParam('xwindow',3)
    iraf.imfilter.median.setParam('ywindow',3)
    iraf.imfilter.median()
  iraf.set(clobber = 'no')
  
  
  return normalization_factor

def illumination_corrections(images):
  for key in OTA_dictionary:
    correction_image = OTA_dictionary[key]+'.med.fits'
    for i in range(len(images)):
      image_to_correct = images[i]+'['+str(key)+']'
      corrected_image = 'illcor_'+OTA_dictionary[key]+'.'+str(images[i][-27:49])
      #print image_to_correct,correction_image,corrected_image
      iraf.unlearn(iraf.imutil.imarith,iraf.imfilter.median)
      iraf.imutil.imarith.setParam('operand1',image_to_correct)
      iraf.imutil.imarith.setParam('op','/')
      iraf.imutil.imarith.setParam('operand2',correction_image)
      iraf.imutil.imarith.setParam('result',illcorpath+corrected_image)
      iraf.imutil.imarith(mode='h')

#make_bpms(images)
#imcombine_lists(images)
#dark_sky_flat()
#normfactor = illumination_corrections(images)


  