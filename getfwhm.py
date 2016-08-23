import sys, os, glob, string
import numpy as np
import astropy as ast
import matplotlib.pyplot as plt
from pyraf import iraf
from tqdm import tqdm
import odi_config as odi

def getfwhm_ota(img, ota, radius=4.0, buff=7.0, width=5.0):
    """
    Get a fwhm estimate for a single OTA using the SDSS catalog stars and
    IRAF imexam (SLOW, but works). Adapted from Kathy Rohde's getfwhm script
    (this implementation is simpler in practice). The radius, buff, and width
    parameters are for the pyraf task rimexam. This fwhm measure comes from
    a gaussian fittype.

    The positions of the SDSS starts are pulled from a ``coords`` file. This
    module automatically fetches the ``coords`` file for the ``img`` and ``ota``
    being processed from the appropriate directory.

    In addition to a median fwhm measurement this module will also
    produce an ouputfile where the positions and fwhm of each source are stored.
    This ``output`` file is used in other modules in the ``odi-tools`` software.
    The name of this ``output`` file is generated based on the ``img`` and
    ``ota``.

    Parameters
    -----------
    img : str
       String containing name of the image currently in use

    ota : str
       Name of ota extension to be used (e.g. OTA33.SCI)

    Returns
    --------
    gfwhm : float
          Median fwhm measure of sources found in the ota field.

    Examples
    --------
    >>> img = 'img1.fits'
    >>> ota = 'OTA33.SCI'
    >>> gfwhm = getfwhm_ota(img,ota)

    """
    # coords= img[0:-5]+'.'+ota+'.sdssxy'
    image = odi.reprojpath+'reproj_'+ota+'.'+str(img[16:])
    coords = odi.coordspath+'reproj_'+ota+'.'+str(img[16:-5])+'.sdssxy'
    print image, coords
    outputfile = odi.coordspath+img[0:-5]+'.'+ota+'.fwhm.log'

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
    print 'median gwfhm in ota',ota+': ',np.median(gfwhm[np.where(gfwhm < 900.0)]),'pixels'# (determined via QR)'
    return np.median(gfwhm[np.where(gfwhm < 900.0)])

def getfwhm_full(img, radius=4.0, buff=7.0, width=5.0):
    """
    Get a fwhm estimate for a stacked image using the SDSS catalog stars and
    IRAF imexam (SLOW, but works). Adapted from Kathy Rohde's getfwhm script
    (this implementation is simpler in practice). The radius, buff, and width
    parameters are for the pyraf task rimexam. This fwhm measure comes from
    a gaussian fittype.

    The positions of the SDSS starts are pulled from a ``coords`` file. This
    module automatically fetches the ``coords`` file for the ``img`` and ``ota``
    being processed from the appropriate directory.

    In addition to a median fwhm measurement this module will also
    produce an ouputfile where the positions and fwhm of each source are stored.
    This ``output`` file is used in other modules in the ``odi-tools`` software.
    The name of this ``output`` file is generated based on the ``img``.

    Parameters
    -----------
    img : str
       String containing name of the image currently in use


    Returns
    --------
    gfwhm : float
          Median fwhm measure of sources found in the ota field.

    Examples
    --------
    >>> img = 'stack1.fits'
    >>> gfwhm = getfwhm_full(img)

    """
    coords = img[:-5]+'.sdssxy'
    outputfile = img[0:-5]+'.fwhm.log'

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
