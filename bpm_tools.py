import os
import sys
import odi_config as odi
from pyraf import iraf

def force_update_bpm(img, ota): 
	image = odi.scaledpath+'scaled_'+ota+'.'+str(img[16:])
	mask_name = odi.bppath+'reproj_mask_'+ota+'.'+str(img[16:])
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
	return

def make_bpms(img, ota):
	# for i in range(len(images)):
	#   for key in OTA_dictionary:
	# create string for mask fits name
	mask_name = odi.bppath+'mask_'+ota+'.'+str(img[16:])
	BPM = mask_name.replace('fits','pl')
	if not os.path.isfile(BPM):
		mask,gaps = odi.mask_ota(img,ota)
		hdu = odi.fits.PrimaryHDU(mask.astype(float))
	if not os.path.isfile(mask_name):
		hdu.writeto(mask_name,clobber=True)
	if not os.path.isfile(mask_name.replace('fits','pl')):
		iraf.unlearn(iraf.imutil.imcopy)
		iraf.imutil.imcopy.setParam('input',mask_name)
		iraf.imutil.imcopy.setParam('output',mask_name.replace('fits','pl'))
		iraf.imutil.imcopy.setParam('verbose','no')
		iraf.imutil.imcopy(mode='h')
	iraf.unlearn(iraf.imutil.hedit)
	iraf.imutil.hedit.setParam('images',img+'['+ota+']')
	iraf.imutil.hedit.setParam('fields','BPM')
	iraf.imutil.hedit.setParam('value',BPM)
	iraf.imutil.hedit.setParam('add','yes')
	iraf.imutil.hedit.setParam('addonly','no')
	iraf.imutil.hedit.setParam('verify','no')
	iraf.imutil.hedit.setParam('update','yes')
	iraf.imutil.hedit(show='no', mode='h')
	return

def main():
    pass

if __name__ == '__main__':
    main()

