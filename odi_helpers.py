import sys, os, glob, string
import numpy as np
import astropy as ast
import matplotlib.pyplot as plt
from pyraf import iraf

import odi_config as odi

def deg_to_sex(ra, dec):
	from astropy import units as u
	from astropy.coordinates import Angle
	rad = Angle(ra * u.deg)
	decd = Angle(dec * u.deg)

	ra = rad.to_string(unit=u.hour, sep=':')
	dec = decd.to_string(unit=u.deg, sep=':')

	return ra, dec

def get_targ_ra_dec(img, ota):
	from astropy.io import fits
	hdulist = fits.open(img)
	hdu = hdulist[0]
	# hdulist.info()
	# print hdu.header
	ra = hdu.header['RA']
	dec = hdu.header['DEC']

	return ra, dec
	hdulist.close()
	
def imcombine_lists(images, filters):
	from astropy.io import fits
	for filter in filters:
		for key in odi.OTA_dictionary:
			list_name =  open(odi.OTA_dictionary[key]+'.'+filter+'.lis',"w")
			for i in range(len(images)):
				hdulist = fits.open(images[i])
				hdr = hdulist[0].header
				filt = hdr['filter']
				if filt == filter:
					print >>list_name,images[i]+'['+str(key)+']'
				hdulist.close()
			list_name.close()
	return	
	
def reproject_ota(img, ota, rad, decd):
	from pyraf import iraf
	image = odi.illcorpath+'illcor_'+ota+'.'+str(img[16:])
	imout = odi.reprojpath+'reproj_'+ota+'.'+str(img[16:])
	iraf.mscred(_doprint=0)
	iraf.clobber='no'
	iraf.unlearn(iraf.mscred.mscimage)
	iraf.mscred.mscimage.format='image'
	iraf.mscred.mscimage.pixmask='yes'
	iraf.mscred.mscimage.verbose='yes'
	iraf.mscred.mscimage.wcssour='parameters'
	iraf.mscred.mscimage.ref=''
	iraf.mscred.mscimage.ra=rad
	iraf.mscred.mscimage.dec=decd
	iraf.mscred.mscimage.scale=0.11
	iraf.mscred.mscimage.rotation=0.0
	iraf.mscred.mscimage.blank=-999
	iraf.mscred.mscimage.interpo='poly5'
	iraf.mscred.mscimage.minterp='poly5'
	iraf.mscred.mscimage.nxbl=2048
	iraf.mscred.mscimage.nybl=2048
	iraf.mscred.mscimage.fluxcon='yes'
	iraf.mscred.mscimage(image,imout)
	
	return

def bgsub_ota(img, ota, apply=False):
	from pyraf import iraf
	image = odi.reprojpath+'reproj_'+ota+'.'+str(img[16:])
	imout = odi.bgsubpath+'bgsub_'+ota+'.'+str(img[16:])
	bg_mean, bg_median, bg_std = odi.mask_ota(img, ota, reproj=True)
	if apply:
	# print bg_mean, bg_median, bg_std
		iraf.unlearn(iraf.imutil.imarith)
		iraf.imutil.imarith.setParam('operand1',image)
		iraf.imutil.imarith.setParam('op','-')
		iraf.imutil.imarith.setParam('operand2',bg_median)
		iraf.imutil.imarith.setParam('result',imout)
		iraf.imutil.imarith.setParam('verbose','yes')
		iraf.imutil.imarith(mode='h')
	return bg_mean, bg_median, bg_std

def stack_otas(ota):
	from pyraf import iraf
	iraf.immatch.imcombine(odi.scaledpath+'scaled_'+ota+'*.fits', odi.otastackpath+ota+'_stack.fits', combine='average', reject='none', offsets='wcs', masktype='goodvalue', maskval=0, blank=-999, scale='none', zero='none', lthresh=-900, hthresh=60000, logfile=ota+'_stack.log')   

def deep_obj_mask(img, ota, apply=False):
	from astropy.io import fits
	from astropy.stats import sigma_clipped_stats
	image = odi.scaledpath+'scaled_'+ota+'.'+str(img[16:])
	ota_mask = 'objmask_'+ota+'.'+str(img[16:17])+'.fits'
	hdulist = fits.open(image)
	hdu_ota = hdulist[0]
	# maskhdu = fits.open(bppath+ota_mask)
	gapshdu = fits.open(odi.bppath+'reproj_mask_'+ota+'.'+str(img[16:]))
	total_mask = gapshdu[0].data
	#maskhdu[0].data + 

	nx, ny = hdu_ota.data.shape
	print nx, ny
	mean1, median1, std1 = sigma_clipped_stats(hdu_ota.data[0:ny/2,0:nx/2], mask=total_mask[0:ny/2,0:nx/2], sigma=3.0, iters=3)
	mean2, median2, std2 = sigma_clipped_stats(hdu_ota.data[0:ny/2,nx/2:nx], mask=total_mask[0:ny/2,nx/2:nx], sigma=3.0, iters=3)
	mean3, median3, std3 = sigma_clipped_stats(hdu_ota.data[ny/2:ny,0:nx/2], mask=total_mask[ny/2:ny,0:nx/2], sigma=3.0, iters=3)
	mean4, median4, std4 = sigma_clipped_stats(hdu_ota.data[ny/2:ny,nx/2:nx], mask=total_mask[ny/2:ny,nx/2:nx], sigma=3.0, iters=3)
	mean = [mean1, mean2, mean3, mean4]
	median = [median1, median2, median3, median4]
	std = [std1, std2, std3, std4]
	# plt.imshow(hdu_ota.data, origin='lower', cmap='Greys_r', vmin=-10., vmax=500.)
	# plt.imshow(total_mask, cmap=random_cmap(), alpha=0.5)
	# plt.show()
	return mean, median, std

def stack_images(refimg):
	from astropy.io import fits
	from pyraf import iraf
	print refimg
	fitsref = fits.open(refimg)
	hduref = fitsref[0]
	objname = hduref.header['object']
	filter_name = hduref.header['filter']
	sky_med = hduref.header['skybg']
	output = objname+'_'+filter_name+'.fits'
	output_bpm = objname+'_'+filter_name+'_bpm.pl'
	iraf.unlearn(iraf.immatch.imcombine, iraf.imutil.imarith)
	iraf.immatch.imcombine(odi.scaledpath+'*'+filter_name+'*.fits', 'temp', combine='average', reject='none', offsets='wcs', masktype='goodvalue', maskval=0, blank=-999, scale='none', zero='none', lthresh=-900, hthresh=60000)
	# iraf.imutil.imarith.setParam('operand1','temp')
	# iraf.imutil.imarith.setParam('op','+')
	# iraf.imutil.imarith.setParam('operand2',sky_med)
	# iraf.imutil.imarith.setParam('result',output)
	# iraf.imutil.imarith.setParam('verbose','yes')
	# iraf.imutil.imarith(mode='h')
	iraf.imutil.imexpr('(a != -999) ? a + b : -999',output,'temp.fits',sky_med)
	iraf.imutil.imexpr('a < 0',output_bpm, output)
	iraf.imutil.imdelete('temp', verify='no')
	iraf.unlearn(iraf.imutil.hedit)
	iraf.imutil.hedit.setParam('images',output)
	iraf.imutil.hedit.setParam('fields','BPM')
	iraf.imutil.hedit.setParam('value',output_bpm)
	iraf.imutil.hedit.setParam('add','yes')
	iraf.imutil.hedit.setParam('addonly','no')
	iraf.imutil.hedit.setParam('verify','no')
	iraf.imutil.hedit.setParam('update','yes')
	iraf.imutil.hedit(show='no', mode='h')
	# iraf.immatch.imcombine(reprojpath+'*.fits', 'test', expm='exp.pl', combine='average', reject='none', offsets='wcs', masktype='goodvalue', maskval=0, blank=-999, scale='none', zero='none', lthresh=-900, hthresh=60000)  
	return output

def instument(img):
    from astropy.io import fits
    """
    A function to grab what version of odi has been used
    should return 'podi' or '5odi'
    Add a line to this effect before all of the other
    odi_process.py stuff
    inst = odi.instument(images[0])
    """
    hdulist = fits.open(img)
    instument_name = hdulist[0].header['INSTRUME']
    hdulist.close()
    print 'Setting instrument to: ', instument_name 
    return instument_name

def main():
    pass

if __name__ == '__main__':
    main()

