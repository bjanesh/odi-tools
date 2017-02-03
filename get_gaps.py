import sys, os, glob, string
import numpy as np
import astropy as ast
import matplotlib.pyplot as plt
from pyraf import iraf
from tqdm import tqdm
import odi_config as odi

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
            hdu.writeto(mask_name,clobber=True)
    if not os.path.isfile(mask_name.replace('fits','pl')):
        iraf.unlearn(iraf.imutil.imcopy)
        iraf.imutil.imcopy.setParam('input',mask_name)
        iraf.imutil.imcopy.setParam('output',mask_name.replace('fits','pl'))
        iraf.imutil.imcopy.setParam('verbose','no')
        iraf.imutil.imcopy(mode='h')
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


def main():
    pass

if __name__ == '__main__':
    main()
