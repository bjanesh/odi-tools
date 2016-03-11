import sys, os, glob, string
import numpy as np
import astropy as ast
import matplotlib.pyplot as plt
from pyraf import iraf
from tqdm import tqdm
import odi_config as odi

def getfwhm_ota(img, ota, radius=4.0, buff=7.0, width=5.0):
	'''
	Get a fwhm estimate for the image using the SDSS catalog stars and IRAF imexam (SLOW, but works)
	Adapted from Kathy's getfwhm script (this implementation is simpler in practice)
	'''
	# coords= img[0:-5]+'.'+ota+'.sdssxy'
	image = odi.reprojpath+'reproj_'+ota+'.'+str(img[16:])
	coords = odi.coordspath+'reproj_'+ota+'.'+str(img[16:-5])+'.sdssxy'
	print image, coords
	outputfile = odi.coordspath+img[0:-5]+'.'+ota+'.fwhm.log'
	# iraf.noao(_doprint=0)
	# iraf.obsutil(_doprint=0)
	# iraf.unlearn("psfmeasure")
	# iraf.noao.obsutil.psfmeasure.setParam('radius',radius)
	# iraf.noao.obsutil.psfmeasure.setParam('iterations',3)
	# iraf.noao.obsutil.psfmeasure.setParam('sbuffer',buff)
	# iraf.noao.obsutil.psfmeasure.setParam('swidth',width)
	# iraf.noao.obsutil.psfmeasure.setParam('display','no')
	# iraf.noao.obsutil.psfmeasure.setParam('frame',1)
	# 
	# if not os.path.isfile(outputfile):
	#     iraf.noao.obsutil.psfmeasure(img, coords='markall', size='GFWHM', imagecur = coordspath+coords, wcs = "logical", logfile=outputfile, StdoutG='/dev/null')

	iraf.tv.rimexam.setParam('radius',radius)
	iraf.tv.rimexam.setParam('buffer',buff)
	iraf.tv.rimexam.setParam('width',width)
	iraf.tv.rimexam.setParam('rplot',20.)
	iraf.tv.rimexam.setParam('center','yes')
	# fit a gaussian, rather than a moffat profile (it's more robust for faint sources)
	iraf.tv.rimexam.setParam('fittype','gaussian')
	iraf.tv.rimexam.setParam('iterati',1)

	if not os.path.isfile(outputfile):
		iraf.tv.imexamine(image, frame=10, logfile = outputfile, keeplog = 'yes', defkey = "a", nframes=0, imagecur = coords,use_display='no',  StdoutG='/dev/null',mode='h')
	# 
	# # unfortunately we have to toss the first measured fwhm value from the median because of the file format    
	# # gfwhm = np.genfromtxt(outputfile, usecols=(3,), skip_header=4, skip_footer=3, unpack=True)
	outputfile_clean = open(outputfile.replace('.log','_clean.log'),"w")
	for line in open(outputfile,"r"):
	    if not 'INDEF' in line:
		outputfile_clean.write(line)
	    if 'INDEF' in line:
		outputfile_clean.write(line.replace('INDEF','999'))
	outputfile_clean.close()
	os.rename(outputfile.replace('.log','_clean.log'),outputfile)
	gfwhm = np.loadtxt(outputfile, usecols=(10,), unpack=True)
	# hdulist = ast.io.fits.open(image)
	# seeing = hdulist[0].header['FWHMSTAR']
	# gfwhm = seeing/0.11
	print 'median gwfhm in ota',ota+': ',np.median(gfwhm),'pixels'# (determined via QR)'
	return np.median(gfwhm)

def getfwhm_full(img, radius=4.0, buff=7.0, width=5.0):
	'''
	Get a fwhm estimate for the image using the SDSS catalog stars and IRAF imexam (SLOW, but works)
	Adapted from Kathy's getfwhm script (this implementation is simpler in practice)
	'''
	coords = img[:-5]+'.sdssxy'
	# print image, coords
	outputfile = img[0:-5]+'.fwhm.log'
	# iraf.noao(_doprint=0)
	# iraf.obsutil(_doprint=0)
	# iraf.unlearn("psfmeasure")
	# iraf.noao.obsutil.psfmeasure.setParam('radius',radius)
	# iraf.noao.obsutil.psfmeasure.setParam('iterations',3)
	# iraf.noao.obsutil.psfmeasure.setParam('sbuffer',buff)
	# iraf.noao.obsutil.psfmeasure.setParam('swidth',width)
	# iraf.noao.obsutil.psfmeasure.setParam('display','no')
	# iraf.noao.obsutil.psfmeasure.setParam('frame',1)
	# 
	# if not os.path.isfile(outputfile):
	#     iraf.noao.obsutil.psfmeasure(img, coords='markall', size='GFWHM', imagecur = coordspath+coords, wcs = "logical", logfile=outputfile, StdoutG='/dev/null')

	iraf.tv.rimexam.setParam('radius',radius)
	iraf.tv.rimexam.setParam('buffer',buff)
	iraf.tv.rimexam.setParam('width',width)
	iraf.tv.rimexam.setParam('rplot',20.)
	iraf.tv.rimexam.setParam('center','yes')
	# fit a gaussian, rather than a moffat profile (it's more robust for faint sources)
	iraf.tv.rimexam.setParam('fittype','gaussian')
	iraf.tv.rimexam.setParam('iterati',1)

	if not os.path.isfile(outputfile):
		iraf.tv.imexamine(img, frame=10, logfile = outputfile, keeplog = 'yes', defkey = "a", nframes=0, imagecur = coords, wcs = "logical", use_display='no',  StdoutG='/dev/null',mode='h')
	outputfile_clean = open(outputfile.replace('.log','_clean.log'),"w")
	for line in open(outputfile,"r"):
	    if not 'INDEF' in line:
		outputfile_clean.write(line)
	    if 'INDEF' in line:
		outputfile_clean.write(line.replace('INDEF','999'))
	outputfile_clean.close()
	os.rename(outputfile.replace('.log','_clean.log'),outputfile)
	# 
	# # unfortunately we have to toss the first measured fwhm value from the median because of the file format    
	# # gfwhm = np.genfromtxt(outputfile, usecols=(3,), skip_header=4, skip_footer=3, unpack=True)
	gfwhm = np.loadtxt(outputfile, usecols=(10,), unpack=True)
	# hdulist = ast.io.fits.open(image)
	# seeing = hdulist[0].header['FWHMSTAR']
	# gfwhm = seeing/0.11
	print 'median gwfhm in ',img+': ',np.median(gfwhm),'pixels'# (determined via QR)'
	return np.median(gfwhm)

def main():
    pass

if __name__ == '__main__':
    main()

