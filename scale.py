import os
import sys
import numpy as np
from astropy.io import fits
from pyraf import iraf
import warnings
import odi_config as odi

def find_ref_image(images):
    # use photometry over the WHOLE image to calculate the scaling factors
    imgs, fwhm, zp_med, zp_std, bg_mean, bg_median, bg_std = np.loadtxt('derived_props.txt', usecols=(0,3,4,5,6,7,8), unpack=True)
    filter_string = np.loadtxt('derived_props.txt', usecols=(2,), unpack=True,dtype=str)
    
    lvls = []
    ams = []
    zps = []
    #print images
    print '#       bg        airmass       zp       zp_std        n_zps'
    for j,im in enumerate(images):
        hdulist = fits.open(im.f)
        airmass = hdulist[0].header['AIRMASS']
        filter  = hdulist[0].header['FILTER']
        these = np.where((imgs.astype(int)==int(im[16])) & (filter_string == filter))
        zp_filter = np.where((imgs.astype(int)==int(im[16])) & (filter_string == filter) & (zp_med < 900.0))
        bg_lvl = np.mean(bg_median[these])
        warnings.simplefilter("error")
        try:
            zp_lvl = np.mean(zp_med[zp_filter])
            zp_lvl_std = np.std(zp_med[zp_filter])
        except RuntimeWarning:
            zp_lvl = 999
            zp_lvl_std = 999
        lvls.append(bg_lvl)
        ams.append(airmass)
        zps.append(zp_lvl)
        hdulist.close()
        print im[16], '%10.3f'%bg_lvl, '%10.3f'%airmass,'%10.3f'%zp_lvl,'%10.3f'%zp_lvl_std,'%10d'%int(len(zp_filter[0]))
        ref_img = np.argmin(np.array(ams))
    print 'reference image:',images[ref_img]
    print np.argmin(np.array(zps))
    return ref_img

def getscale(images, refimg, verbose=True):
	# img = bgsubpath+'bgsub_'+ota+'.'+img.stem()
	# refimg = bgsubpath+'bgsub_'+ota+'.'+str(refref_img.stem())
	scale = {}
	std = {}
	sigThreshold = 0.006
	for img in images:
		print 'Calculating scaling factor for', img
		fitsref = fits.open(refimg)
		hduref = fitsref[0]
		fitsimg = fits.open(img.f)
		hduimg = fitsimg[0]

		exp_ref = hduref.header['EXPTIME']
		exp_img = hduimg.header['EXPTIME']

		n = 1
		expRatio = (float(exp_ref)/float(exp_img))

		# use photometry over the whole image, so read in all the phot output files and append

		for key in odi.OTA_dictionary.keys():
			ota = odi.OTA_dictionary[key] 
			img_phot = odi.coordspath+img.nofits()+'.'+ota+'.sdssphot'
			ref_phot = odi.coordspath+refimg.nofits()+'.'+ota+'.sdssphot'
			img_peak = odi.coordspath+img.nofits()+'.'+ota+'.fwhm.log'
			ref_peak = odi.coordspath+refimg.nofits()+'.'+ota+'.fwhm.log'
			img_coords = odi.coordspath+'reproj_'+ota+'.'+img.base()+'.sdssxy'
			ref_coords = odi.coordspath+'reproj_'+ota+'.'+str(refrefimg.base())+'.sdssxy'
			# print 'getting data from'
			# print img_phot
			# print ref_phot
			print ota, img
			if key==1:
				# read in the instrumental magnitudes from the appropriate phot output files
				magTempA, magErrTempA, skyTempA = np.loadtxt(img_phot, usecols=(1,2,3), unpack=True)
				idA = np.loadtxt(img_phot, usecols=(0,), unpack=True)
				magTempRef, magErrTempRef, skyTempRef = np.loadtxt(ref_phot, usecols=(1,2,3), unpack=True)
				idRef = np.loadtxt(ref_phot, usecols=(0,), unpack=True)
				peakTempA = np.loadtxt(img_peak, usecols=(9,), unpack=True)
				peakTempRef = np.loadtxt(ref_peak, usecols=(9,), unpack=True)
				raTempA, decTempA = np.loadtxt(img_coords, usecols=(2,3), unpack=True)
				raTempRef, decTempRef = np.loadtxt(ref_coords, usecols=(2,3), unpack=True)
				uniqidA = zip(raTempA,decTempA)
				uniqidRef = zip(raTempRef,decTempRef)

				# gID_keep = gID - 1
				# iID_keep = iID - 1
				# keep = list(set(gID_keep).intersection(iID_keep))

				# keep the common elements between the two phot lists using their list index
				keepRef = [i for i,element in enumerate(uniqidRef) if element in uniqidA]
				keepA = [i for i,element in enumerate(uniqidA) if element in uniqidRef]
				# print keepRef, keepA
				# kept = list(set(keptRef).intersection(keptA))
				# kept = list(np.array(kept).astype(int))
				# keepRef = [i for i,element in enumerate(idRef) if element in kept]
				# keepA = [i for i,element in enumerate(idA) if element in kept]
				# print kept, keepRef, keepA

				magA = magTempA[np.array(keepA)]
				magRef = magTempRef[np.array(keepRef)]
				magErrA = magErrTempA[np.array(keepA)]
				magErrRef = magErrTempRef[np.array(keepRef)]
				skyA = skyTempA[np.array(keepA)]
				skyRef = skyTempRef[np.array(keepRef)]
				peakA = peakTempA[np.array(keepA)]
				peakRef = peakTempRef[np.array(keepRef)]
				# for k,m in enumerate(magA):
				#     print peakA[k], peakRef[k], magA[k], magRef[k]

				# and then get rid of poorly measuerd objects probably in cell gaps
				keep = np.where((skyA>0.0) & (skyRef > 0.0) & (10000.0<peakA) & (10000.0<peakRef) & (45000.0>peakA) & (45000.0>peakRef))
				# (magErrA<0.01) & (magErrRef<0.01) &
				magA = magA[keep]
				magRef = magRef[keep]
			else:
				# read in the instrumental magnitudes from the appropriate phot output files
				magTempB, magErrTempB, skyTempB = np.loadtxt(img_phot, usecols=(1,2,3), unpack=True)
				idB = np.loadtxt(img_phot, usecols=(0,), unpack=True)
				magTempRefB, magErrTempRefB, skyTempRefB = np.loadtxt(ref_phot, usecols=(1,2,3), unpack=True)
				idRefB = np.loadtxt(ref_phot, usecols=(0,), unpack=True)
				peakTempB = np.loadtxt(img_peak, usecols=(9,), unpack=True)
				peakTempRefB = np.loadtxt(ref_peak, usecols=(9,), unpack=True)
				raTempB, decTempB = np.loadtxt(img_coords, usecols=(2,3), unpack=True)
				raTempRefB, decTempRefB = np.loadtxt(ref_coords, usecols=(2,3), unpack=True)
				uniqidB = zip(raTempB,decTempB)
				uniqidRefB = zip(raTempRefB,decTempRefB)

				# gID_keep = gID - 1
				# iID_keep = iID - 1
				# keep = list(set(gID_keep).intersection(iID_keep))

				# keep the common elements between the two phot lists using their list index
				keepRefB = [i for i,element in enumerate(uniqidRefB) if element in uniqidB]
				keepB = [i for i,element in enumerate(uniqidB) if element in uniqidRefB]
				# print keepRef, keepA
				# kept = list(set(keptRef).intersection(keptA))
				# kept = list(np.array(kept).astype(int))
				# keepRef = [i for i,element in enumerate(idRef) if element in kept]
				# keepA = [i for i,element in enumerate(idA) if element in kept]
				# print kept, keepRef, keepA

				magB = magTempB[np.array(keepB)]
				magRefB = magTempRefB[np.array(keepRefB)]
				magErrB = magErrTempB[np.array(keepB)]
				magErrRefB = magErrTempRefB[np.array(keepRefB)]
				skyB = skyTempB[np.array(keepB)]
				skyRefB = skyTempRefB[np.array(keepRefB)]
				peakB = peakTempB[np.array(keepB)]
				peakRefB = peakTempRefB[np.array(keepRefB)]

				# for k,m in enumerate(magB):
				#     print peakB[k], peakRefB[k], magB[k], magRefB[k]
				# and then get rid of poorly measuerd objects probably in cell gaps
				keep = np.where((skyB>0.0) & (skyRefB > 0.0) & (10000.0<peakB) & (10000.0<peakRefB) & (45000.0>peakB) & (45000.0>peakRefB))
				# (magErrB<0.01) & (magErrRefB<0.01) & 
				magB = magB[keep]
				magRefB = magRefB[keep]
				magA = np.hstack((magA, magB))
				magRef = np.hstack((magRef, magRefB))

		rat = np.power(10.0,-0.4*(magA-magRef))/expRatio
		if verbose:
			for i,r in enumerate(rat):
				print img, i, r

		sigTest = np.std(rat)
		print len(rat), np.mean(rat), np.median(rat), np.std(rat), n
		if sigTest <= sigThreshold:
			scale[img] = np.mean(rat)
			std[img] = np.std(rat)
		else:
			while sigTest > sigThreshold:
				magTempA = magA
				magTempRef = magRef
				magA = magTempA[np.where(abs(rat-np.median(rat))<sigTest)]
				magRef = magTempRef[np.where(abs(rat-np.median(rat))<sigTest)]
				rat = np.power(10.0,-0.4*(magA-magRef))/expRatio
				for i,r in enumerate(rat):
					print magA[i], magRef[i], r
				sigTest = np.std(rat)

				n = n + 1
				if n > 10:
					print "Iteration did not converge to sigma <", repr(sigThreshold),"for", img
					print "Quitting..."
					exit()
				print len(rat), np.mean(rat), np.median(rat), np.std(rat), n
			scale[img] = np.mean(rat)
			std[img] = np.std(rat)

	return scale, std

def scale_ota(img, ota, scale):
	image = odi.bgsubpath+'bgsub_'+ota+'.'+img.stem()
	imout = odi.scaledpath+'scaled_'+ota+'.'+img.stem()
	iraf.unlearn(iraf.imutil.imarith)
	iraf.imutil.imarith.setParam('operand1',image)
	iraf.imutil.imarith.setParam('op','/')
	iraf.imutil.imarith.setParam('operand2',scale)
	iraf.imutil.imarith.setParam('result',imout)
	iraf.imutil.imarith.setParam('verbose','yes')
	iraf.imutil.imarith(mode='h')
	return 


def main():
    pass

if __name__ == '__main__':
    main()

