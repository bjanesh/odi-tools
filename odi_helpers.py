import sys, os, glob, string
import numpy as np
import astropy as ast
import matplotlib.pyplot as plt
from pyraf import iraf
from tqdm import tqdm

import odi_config as odi

def deg_to_sex(ra, dec):
    """
    Convert an Ra and Dec position in decimal degrees to hexadecimal
    Parameters
    ----------
    ra : float
        Ra in decimal degrees
    dec : float
        Dec in decimal degrees

    Returns
    -------
    ra : str
        Ra in hexadecimal HH:MM:SS

    dec : str
        Dec in hexadecimal DD:MM:SS
    """
    from astropy import units as u
    from astropy.coordinates import Angle
    rad = Angle(ra * u.deg)
    decd = Angle(dec * u.deg)

    ra = rad.to_string(unit=u.hour, sep=':')
    dec = decd.to_string(unit=u.deg, sep=':')

    return ra, dec

def get_targ_ra_dec(img, ota):
    """
    Get and return the Ra and Dec of an OTA based on the ``RA`` and ``DEC``
    header cards.

    Parameters
    ----------
    img : str
        Name of image
    ota : str
        Name of ota

    Returns
    -------
    ra : float
        Ra in decimal degrees

    dec : str
        Dec in decimal degrees
    """
    from astropy.io import fits
    hdulist = fits.open(img.f)
    hdu = hdulist[0]
    # hdulist.info()
    # print hdu.header
    ra = hdu.header['RA']
    dec = hdu.header['DEC']

    return ra, dec
    hdulist.close()
    
def list_wcs_coords(img, ota, gapmask, inst,output='radec.coo', gmaglim=20., stars_only=True, offline = False, source = 'sdss'):
    """
    Create the files needed to fix the WCS solution on a given ota. This
    function will create lists of SDSS, 2MASS, or Gaia sources depending on the
    options selected by the user. If this function is run in the ``offline``
    mode, the source catalogs will be taken from the files produced by
    :py:func:`offlinecats`

    Parameters
    ----------
    img : str
        Name of image
    ota : str
        Name of OTA
    gapmask : numpy array
        A numpy array of the gap location on the ota. This can be produced by
        the function :py:func:`get_gaps.get_gaps`. The gap mask is used to
        filter out stars that fall in or near gaps on the ota.
    int : str
        Version of ODI used, ``podi`` or ``5odi``
    output : str
        Desired name of the output catalog
    gmaglim : float
        Magnitude limit in the g band for sources that will be included in the
        catalogs. This might need to be adjusted according to your data. If it
        is too bright, there might not be enough sources to produces a good
        WCS solution.
    stars_only : bool
        When using SDSS sources this only includes sources flagged as stars
    offline : bool
        When ``True`` this function will use the catalogs produced by
        :py:func:`offlinecats`. If ``False`` this function will query the
        online ``SDSS`` catalog for sources.
    source : str
        Name of desired source catalog. Must be either ``sdss``,``twomass``, or
        ``gaia``.

    Note
    ----
    This functions produces three files for each ota in the ``coords``
    directory with the following naming scheme:

    1. ``img.nofits()+'.'+ota+'.radec.coo'``
    2. ``img.nofits()+'.'+ota+'.radec.coo.px'``
    3. ``img.nofits()+'.'+ota+'.sdssxy'``
    """

    if offline == False:
        xdim, ydim = odi.get_sdss_coords(img, ota, inst,output=odi.coordspath+img.nofits()+'.'+ota+'.sdss')
        ras,decs,psfMag_u,psfMagErr_u,psfMag_g,psfMagErr_g,psfMag_r,psfMagErr_r,psfMag_i,psfMagErr_i,psfMag_z,psfMagErr_z = np.loadtxt(odi.coordspath+img.nofits()+'.'+ota+'.sdss',usecols=(0,1,2,3,4,5,6,7,8,9,10,11), unpack=True, delimiter=',', skiprows=2)
        probPSF = np.loadtxt(odi.coordspath+img.nofits()+'.'+ota+'.sdss', usecols=(12,), dtype=int, unpack=True, delimiter=',', skiprows=2)
        coords2 = zip(ras[np.where((psfMag_g<gmaglim) & (probPSF==1))],decs[np.where((psfMag_g<gmaglim) & (probPSF==1))])
    if offline == True and source == 'sdss':
        sdss_cat = odi.sdsspath+'offline_'+ota+'.'+img.base()+'.sdss'
        print 'Using Ra and Dec from:', sdss_cat,'for fixwcs'
        ras,decs,psfMag_u,psfMagErr_u,psfMag_g,psfMagErr_g,psfMag_r,psfMagErr_r,psfMag_i,psfMagErr_i,psfMag_z,psfMagErr_z = np.loadtxt(sdss_cat,usecols=(0,1,2,3,4,5,6,7,8,9,10,11), unpack=True, delimiter=',', skiprows=1)
        coords2 = zip(ras[np.where(psfMag_g<gmaglim)],decs[np.where(psfMag_g<gmaglim)])
    if offline == True and source == 'twomass':
        twomass_cat = odi.twomasspath+'offline_'+ota+'.'+img.base()+'.mass'
        ras,decs = np.loadtxt(twomass_cat,usecols=(2,3), unpack=True, delimiter=',', skiprows=1)
        # Just creating dummy variables so that the file formats remain the same for other functions
        psfMag_u       = np.ones(len(ras))
        psfMagErr_u    = np.ones(len(ras))
        psfMag_g       = np.ones(len(ras))
        psfMagErr_g    = np.ones(len(ras))
        psfMag_r       = np.ones(len(ras))
        psfMagErr_r    = np.ones(len(ras))
        psfMag_i       = np.ones(len(ras))
        psfMagErr_i    = np.ones(len(ras))
        psfMag_z       = np.ones(len(ras))
        psfMagErr_z    = np.ones(len(ras))
        coords2 = zip(ras,decs)
    if source == 'gaia':
        gaia_cat = odi.gaiapath+'offline_'+ota+'.'+img.base()+'.gaia'
        ras,decs = np.loadtxt(gaia_cat,usecols=(0,1), unpack=True, delimiter=',', skiprows=1)
        # Just creating dummy variables so that the file formats remain the same
        # for other functions
        psfMag_u       = np.ones(len(ras))
        psfMagErr_u    = np.ones(len(ras))
        psfMag_g       = np.ones(len(ras))
        psfMagErr_g    = np.ones(len(ras))
        psfMag_r       = np.ones(len(ras))
        psfMagErr_r    = np.ones(len(ras))
        psfMag_i       = np.ones(len(ras))
        psfMagErr_i    = np.ones(len(ras))
        psfMag_z       = np.ones(len(ras))
        psfMagErr_z    = np.ones(len(ras))
        coords2 = zip(ras,decs)

    hdulist = odi.fits.open(img.f)
    hdu = odi.tan_header_fix(hdulist[ota])

    if offline == True:
        xdim = hdu.header['NAXIS1']
        ydim = hdu.header['NAXIS2']

    w = odi.WCS(hdu.header)
    pixcrd2 = w.wcs_world2pix(coords2, 1)
    pixid = []
    with open(odi.coordspath+output,'w+') as f:
        with open(odi.coordspath+output+'.pix', 'w+') as fp:
            with open(odi.coordspath+img.nofits()+'.'+ota+'.sdssxy', 'w+') as fxy:
                for i,c in enumerate(coords2):
                    if  20.0 <= pixcrd2[i,0] < xdim-100.0 and 20.0 <= pixcrd2[i,1] < ydim-100.0:
                        # make an image cutout of the gap mask
                        x, y = int(round(pixcrd2[i,0])), int(round(pixcrd2[i,1]))
                        cutout = gapmask[y-30:y+30,x-30:x+30]
                        # print cutout
                        if not (cutout.astype(bool)).any():
                            pixid.append(i)
                            r, d = odi.deg_to_sex(c[0], c[1])
                            print >> f, r, d, psfMag_g[i]
                            print >> fp, pixcrd2[i,0], pixcrd2[i,1], i, 'm'
                            print >> fxy, pixcrd2[i,0], pixcrd2[i,1], ras[i],decs[i],psfMag_u[i],psfMagErr_u[i],psfMag_g[i],psfMagErr_g[i],psfMag_r[i],psfMagErr_r[i],psfMag_i[i],psfMagErr_i[i],psfMag_z[i],psfMagErr_z[i]

    pixid = np.array(pixid)
    pixcrd3 = pixcrd2[pixid]
    hdulist.close()
    return pixcrd3

def fix_wcs(img, ota, coords='radec.coo', iters=3):
    """
    Try to improve the WCS solution of an OTA using the IRAF task ``msccmatch``.
    This function will use the ``img.nofits()+'.'+ota+'.radec.coo'`` file produced
    by :py:func:`list_wcs_coords`.

    Parameters
    ----------
    img : str
        Name of image
    ota : str
        Name of OTA
    coords : str
        Name of coordinate file
    iter : int
        Number of desired iterations of  ``msccmatch``. It is still being
        tested, but one might be all that is necessary, especially if using the
        Gaia catalog.

    Note
    ----
    This function is set up to use the files in the ``illcor`` directory. The
    following are the parameters used by ``msccmatch``.

    - verbose='yes'
    - usebpm='no'
    - nsearch=250
    - search=30
    - rsearch=0.2
    - cfrac=.5
    - csig=0.1
    - nfit=5
    - rms=1.0
    - maxshif=5.0
    - fitgeom="general"
    - update='yes'
    - interac='yes'
    - fit='no',
    - accept='yes'
    - Stdout=1
    """
    image = odi.illcorpath+'illcor_'+ota+'.'+img.stem()
    iraf.mscred(_doprint=0)
    iraf.unlearn(iraf.mscred.msccmatch)

    for i in range(iters):
        fix = iraf.msccmatch(input=image,
                             coords=odi.coordspath+coords,
                             usebpm='no',
                             verbose='yes',
                             nsearch=250,
                             search=30,
                             rsearch=0.2,
                             cfrac=.5,
                             csig=0.1,
                             nfit=5,
                             rms=1.0,
                             maxshif=5.0,
                             fitgeom="general",
                             update='yes',
                             interac='yes',
                             fit='no',
                             accept='yes',
                             Stdout=1)
        print 'fixing WCS for',img.f+'['+ota+'], iter ='+repr(i)
        print fix[-6]
        print fix[-5]
        print fix[-4]
        print fix[-3]
        print fix[-2]

def fix_wcs_full(img, coords='radec.coo', iters=1):
    """
    Try to improve the WCS solution of a final stacked image.

    Parameters
    ----------
    img : str
        Name of image
    coords : str
        Name of coordinate file
    iter : int
        Number of desired iterations of  ``msccmatch``. It is still being
        tested, but one might be all that is necessary, especially if using the
        Gaia catalog.

    Note
    ----
    This function is set up to use the files in the ``illcor`` directory. The
    following are the parameters used by ``msccmatch``.

    - verbose='yes'
    - usebpm='no'
    - nsearch=250
    - search=30
    - rsearch=0.2
    - cfrac=.5
    - csig=0.1
    - nfit=5
    - rms=1.0
    - maxshif=5.0
    - fitgeom="general"
    - update='yes'
    - interac='yes'
    - fit='no',
    - accept='yes'
    - Stdout=1
    """
    print coords
    iraf.mscred(_doprint=0)
    iraf.unlearn(iraf.mscred.msccmatch)
    # otaext = {'33':'[1]','34':'[2]','44':'[3]','43':'[4]','42':'[5]','32':'[6]','22':'[7]','23':'[8]','24':'[9]'}
    for i in range(iters):
        fix = iraf.msccmatch(input=img,
                             coords=coords,
                             usebpm='no',
                             verbose='yes',
                             nsearch=250,
                             search=30,
                             rsearch=0.2,
                             cfrac=.5,
                             csig=0.1,
                             nfit=5,
                             rms=1.0,
                             maxshif=5.0,
                             fitgeom="general",
                             update='yes',
                             interac='yes',
                             fit='no',
                             accept='yes',
                             Stdout=1)
        print 'fixing WCS for',img.f+', iter ='+repr(i)
        print fix[-6]
        print fix[-5]
        print fix[-4]
        print fix[-3]
        print fix[-2]

def repair_bad_wcs(img, ota, refimg, refota):
    print 'repairing bad wcs solution for',img.f+'['+ota+']...'
    # get good CD matrix values from the reference image
    # refimg = refimg.f+'['+refota+']'
    refhdu = odi.fits.open(refimg)
    pvlist = refhdu[refota].header['PV*']
    for pv in pvlist:
        tpv = 'T'+pv
        refhdu[refota].header.rename_keyword(pv, tpv, force=False)
    w_ref = odi.WCS(refhdu[refota].header)
    print w_ref.wcs.cd, w_ref.wcs.crpix, w_ref.wcs.crval

    # get the bad WCS info so we can do some checking
    # image = img.f+'['+ota+']'
    hdu = odi.fits.open(img.f)
    pvlist = hdu[refota].header['PV*']
    for pv in pvlist:
        tpv = 'T'+pv
        hdu[refota].header.rename_keyword(pv, tpv, force=False)
    w = odi.WCS(hdu[ota].header)
    print w.wcs.cd, w.wcs.crpix, w.wcs.crval

def repair_wcs_keywords(img):
    hdulist = odi.fits.open(img.f, mode='update')
    existing_radesys = hdulist[0].header['RADESYS']
    print img.f
    print '--> Existing RADESYS value:', existing_radesys
    correct_radesys = existing_radesys.strip("'").strip()
    print '--> Correct RADESYS value:', correct_radesys
    hdulist[0].header["RADESYS"] = correct_radesys
    # print 'fixing CTYPES in OTA headers'
    for k in tqdm(img.otas):
        existing_ctype1 = hdulist[k].header['CTYPE1']
        existing_ctype2 = hdulist[k].header['CTYPE2']
        correct_ctype1 = existing_ctype1.replace('TAN','TPV')
        correct_ctype2 = existing_ctype2.replace('TAN','TPV')
        hdulist[k].header['CTYPE1'] = correct_ctype1
        hdulist[k].header['CTYPE2'] = correct_ctype2

    hdulist.flush()
    hdulist.close()
    
def getfwhm_ota(img, ota, gaia=False, radius=4.0, buff=7.0, width=5.0):
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
    # coords= img.nofits()+'.'+ota+'.sdssxy'
    image = odi.reprojpath+'reproj_'+ota+'.'+img.stem()
    if gaia:
        coords = odi.coordspath+'reproj_'+ota+'.'+img.base()+'.gaiaxy'
    else:
        coords = odi.coordspath+'reproj_'+ota+'.'+img.base()+'.sdssxy'
    print image, coords
    outputfile = odi.coordspath+img.nofits()+'.'+ota+'.fwhm.log'

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
    coords = img.nofits()+'.sdssxy'
    outputfile = img.nofits()+'.fwhm.log'

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
    print 'median gwfhm in ',img.f+': ',np.median(gfwhm),'pixels'# (determined via QR)'
    return np.median(gfwhm)

def imcombine_lists(images, filters):
    """
    Create a list files for OTAs, sorted by filter, that will be later
    combined to make a dark sky flat.

    Parameters
    ----------
    images : list
        List of images
    filters : list
        List of filters present in the images list

    Note
    ----
    Nothing is returned by this function. A file will be created for each OTA
    following this naming scheme: ``OTA##.SCI.filter.lis``. For example:
    ``OTA22.SCI.odi_g.lis``. Within each of these files will be a list of
    images to combine.
    """
    from astropy.io import fits
    for filter in filters:
        for key in odi.OTA_dictionary:
            list_name =  open(odi.OTA_dictionary[key]+'.'+filter+'.lis',"w")
            for i in range(len(images)):
                hdulist = fits.open(images[i].f)
                hdr = hdulist[0].header
                filt = hdr['filter']
                if filt == filter:
                    print >>list_name,images[i].f+'['+str(key)+']'
                hdulist.close()
            list_name.close()
    return

def reproject_ota(img, ota, rad, decd, wcsref):
    """
    Use the IRAF task ``mscimage`` in the ``mscred`` package to reproject
    an OTA to a reference tangent plane with constant pixel scale.

    Parameters
    ----------
    img : str
        Name of image being processed
    ota : str
        Name of current ``ota`` being processed in ``img``
    rad : float
        Reference Ra position for reprojection
    decd : float
        Reference Ra position for reprojection
    wcfreg : str
        Name of image and ota to be used as the reference image for ``mscimage``

    Note
    ----
    Nothing is returned by this function but the reprojected ota is saved to the
    ``repreopath``. The pipeline is setup to use OTA33 of
    the first image in the images list as the reference image for this function.

    Here is how the ``mscimage`` IRAF parameters are set:

    - iraf.mscred.mscimage.format='image'
    - iraf.mscred.mscimage.pixmask='yes'
    - iraf.mscred.mscimage.verbose='yes'
    - iraf.mscred.mscimage.wcssour='image'
    - iraf.mscred.mscimage.ref=wcsref
    - iraf.mscred.mscimage.ra=rad
    - iraf.mscred.mscimage.dec=decd
    - iraf.mscred.mscimage.scale=0.11
    - iraf.mscred.mscimage.rotation=0.0
    - iraf.mscred.mscimage.blank=-999
    - iraf.mscred.mscimage.interpo='poly5'
    - iraf.mscred.mscimage.minterp='poly5'
    - iraf.mscred.mscimage.nxbl=4096
    - iraf.mscred.mscimage.nybl=4096
    - iraf.mscred.mscimage.fluxcon='yes'
    - iraf.mscred.mscimage(image,imout)

    """
    from pyraf import iraf
    image = odi.illcorpath+'illcor_'+ota+'.'+img.stem()
    imout = odi.reprojpath+'reproj_'+ota+'.'+img.stem()
    iraf.mscred(_doprint=0)
    iraf.clobber='no'
    iraf.unlearn(iraf.mscred.mscimage)
    iraf.mscred.mscimage.format='image'
    iraf.mscred.mscimage.pixmask='yes'
    iraf.mscred.mscimage.verbose='yes'
    iraf.mscred.mscimage.wcssour='image'
    iraf.mscred.mscimage.ref=wcsref
    iraf.mscred.mscimage.ra=rad
    iraf.mscred.mscimage.dec=decd
    iraf.mscred.mscimage.scale=0.11
    iraf.mscred.mscimage.rotation=0.0
    iraf.mscred.mscimage.blank=-999
    iraf.mscred.mscimage.interpo='poly5'
    iraf.mscred.mscimage.minterp='poly5'
    iraf.mscred.mscimage.nxbl=4096
    iraf.mscred.mscimage.nybl=4096
    iraf.mscred.mscimage.fluxcon='yes'
    iraf.mscred.mscimage(image,imout)

    return

def bkg_boxes(hdu,nboxes,length,sources):
    """
    Function to calculate the sigma clipped statistics of a number of randomly
    generated boxes over an ota.

    Parameters
    ----------
    hdu : fits object
        Hdulist that as been opened by astropy.fits.io

    nboxes : int
        Number of random boxes to generate over the ota

    length : int
        Length of side of box in pixels

    sources : bool
        If ``True`` any sources detected in a given box will be masked before
        calculating the background statistics

    Returns
    -------
    bg_stats : numpy array
        Array containing the background stats of all of the boxes
    bg_median : float
        Median background level of the boxes
    med_std : float
        Median standard deviation of the background level in each box
    std_std : float
        Standard deviation of the standard deviations in each box
    centers : list
        Pixel centers of each box

    """
    image = hdu.data

    #Get length of image in each axis
    naxis1 = hdu.header['NAXIS1']
    naxis2 = hdu.header['NAXIS2']

    #generate the centers of n random boxes.

    box_centers = np.random.randint(length,np.min([naxis1-length,naxis2-length]),size=(nboxes,2))

    #divide length by 2
    # another numpy error on wopr (can't convert to integer), so don't worry about integer arithmetic
    side = length/2

    bg_stats = []
    centers = []
    for center in range(len(box_centers)):
        x1 = int(box_centers[center][0]-side)
        x2 = int(box_centers[center][0]+side)
        y1 = int(box_centers[center][1]-side)
        y2 = int(box_centers[center][1]+side)

        #Check to ensure that box is within image
        if (x1 > side and x2 < naxis1-1.5*side) and (y1 > side and y2 < naxis2-1.5*side):
            centers.append(box_centers[center])
            """
            The centers that are within the image bounds are returned
            in case you need to examine the regions used.
            """
            box = image[x1:x2,y1:y2]

            if np.isnan(box).any() == False and (box >= 0).all() == True:
                """
                Only boxes with non-negative values are kept.
                This should help deal with cell gaps
                The sigma and iter values might need some tuning.
                """
                mean, median, std = sigma_clipped_stats(box, sigma=3.0)
                if std >= 2.0*np.sqrt(median):
                    pass
                else:
                    if sources == False:
                        bg_stats.append((mean, median, std))
                    if sources == True:
                        threshold = median + (std * 2.)
                        segm_img = detect_sources(box, threshold, npixels=20)
                        mask = segm_img.data.astype(np.bool)# turn segm_img into a mask
                        selem = np.ones((10, 10))    # dilate using a 25x25 box
                        mask2 = binary_dilation(mask, selem)
                        #new_mask = mask_first_pass + mask2
                        new_mask = mask2

                        mean_mask, median_mask, std_mask = sigma_clipped_stats(box, sigma=3.0, mask=new_mask)
                        bg_stats.append((mean_mask, median_mask, std_mask))

    bg_stats = np.reshape(np.array(bg_stats),(len(bg_stats),3))
    centers = np.reshape(np.array(centers),(len(centers),2))
    #Calculate median std of Background
    med_std = np.median(bg_stats[:,2])
    #calculate standard deviation of the std values
    std_std = np.std(bg_stats[:,2])
    #median
    bg_median = np.median(bg_stats[:,1])

    #Locate the box that had the largest std
    #Array will be returned for plotting if wanted
    # max_std = np.argmax(bg_stats[:,2])
    # max_center = centers[max_std]
    # max_box = image[max_center[0]-side:max_center[0]+side,max_center[1]-side:max_center[1]+side]
    # max_box is currently not needed.
    max_box = 1.0

    return bg_stats,bg_median,med_std,std_std,centers,max_box


def bgsub_ota(img, ota, apply=False):
    """
    Subtract a background level from an OTA.

    Parameters
    ----------
    img : str
        Name of image being processed
    ota : str
        Name of current ``ota`` being processed in ``img``
    apply : bool
        If ``True`` the background level is calculated and subtracted. If
        ``False``, the background level is calculated by not subtracted from
        the current ``ota``.
    Returns
    -------
    bg_mean : float
        Mean background level
    bg_median : float
        Meidan background level
    bg_std : float
        Standard deviation of background

    Note
    ----
    This function calls ``odi.mask_ota`` to calculate the background statistics
    """
    from pyraf import iraf
    image = odi.reprojpath+'reproj_'+ota+'.'+img.stem()
    imout = odi.bgsubpath+'bgsub_'+ota+'.'+img.stem()
    bg_mean, bg_median, bg_std = odi.mask_ota(img, ota, reproj=True)
    tqdm.write('subtracting {:7.2f} from {:s}'.format(bg_median, image))
    if apply:
        # print bg_mean, bg_median, bg_std
        iraf.unlearn(iraf.imutil.imarith)
        iraf.imutil.imarith.setParam('operand1',image)
        iraf.imutil.imarith.setParam('op','-')
        iraf.imutil.imarith.setParam('operand2',bg_median)
        iraf.imutil.imarith.setParam('result',imout)
        iraf.imutil.imarith.setParam('verbose','no')
        iraf.imutil.imarith(mode='h')
    return bg_mean, bg_median, bg_std

def stack_otas(ota):
    """
    (Currently not used)
    Stack the OTAs at the end of the ``odi_process.py`` using the ``IRAF`` total
    ``imcombine``. Here are the the settings used by ``imcombine``:

    - combine='average'
    - reject='none'
    - offsets='wcs'
    - masktype='goodvalue'
    - maskval=0
    - blank=-999
    - scale='none'
    - zero='none'
    - lthresh=-900
    - hthresh=60000
    - logfile=ota+'_stack.log'

    Parameters
    ----------
    ota : str
        Name of current ``ota`` being processed
    """
    from pyraf import iraf
    iraf.immatch.imcombine(odi.scaledpath+'scaled_'+ota+'*.fits', odi.otastackpath+ota+'_stack.fits', combine='average', reject='none', offsets='wcs', masktype='goodvalue', maskval=0, blank=-999, scale='none', zero='none', lthresh=-900, hthresh=60000, logfile=ota+'_stack.log')

def deep_obj_mask(img, ota, apply=False):
    """
    Currently not used
    """
    from astropy.io import fits
    from astropy.stats import sigma_clipped_stats
    image = odi.scaledpath+'scaled_'+ota+'.'+img.stem()
    ota_mask = 'objmask_'+ota+'.'+str(img.dither())+'.fits'
    hdulist = fits.open(image)
    hdu_ota = hdulist[0]
    # maskhdu = fits.open(bppath+ota_mask)
    gapshdu = fits.open(odi.bppath+'reproj_mask_'+ota+'.'+img.stem())
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

def find_new_bg(refimg, filter):
    """
    Calculate a new background level to be added to the OTAs before
    stacking

    Parameters
    ----------
    refimg : str
        Reference image for background calculation
    filter : str
        Filter of the reference image

    Returns
    -------
    sky_med : float
        Median background level

    """
    img, ota, filt, fwhm, zp_med, zp_std, bg_mean, bg_med, bg_std = np.loadtxt('derived_props.txt', usecols=(0,1,2,3,4,5,6,7,8), dtype=str, unpack=True)
    keep = np.where((img == refimg.dither()) & (filt==filter))

    sky_med = np.median(bg_med[keep].astype(float))
    sky_mean = np.median(bg_mean[keep].astype(float))
    sky_std = np.median(bg_std[keep].astype(float))
    print 'calculated sky median, mean, std to re-add:', sky_med, sky_mean, sky_std
    return sky_med, sky_mean, sky_std
    
def is_guide_ota(img, ota):
    """
    Determines whether the specified image OTA was used for guiding.
    
    Parameters
    ----------
    img : str
        Name of image being processed
    ota : str
        Name of current ``ota`` being processed in ``img``

    Returns
    -------
    guide : boolean
        True if guide OTA, False if not
    """
    from astropy.io import fits
    from photutils.segmentation import detect_sources, source_properties, properties_table
    guide = False
    check_ota = 'illcor_'+ota+'.'+img.stem()
    hdu = fits.open(check_ota)
    data = hdu[0].data
    segm = detect_sources(data, 1.0, 50000)
    props = source_properties(data, segm)
    corners = []
    for p in props:
        cutout = p.data_cutout
        yc, xc = int(p.cutout_centroid[0].value), int(p.cutout_centroid[1].value)
        bgcent = cutout[yc,xc]
        bgcorn = np.array([cutout[0,0], cutout[0,-1], cutout[-1,0], cutout[-1,-1]])
        bgrat = bgcorn/bgcent
        corners.append(np.median(bgrat))
    med_ratio = np.median(corners)
    print med_ratio
    if med_ratio > 3.:
        guide = True
    return guide

def make_stack_list(object, filter, inst):
    """
    Makes a list of images to be stacked using ``stack_images()``. This list
    does not include the guiding OTAs as determined by ``derived_props.txt``.

    Parameters
    ----------
    object : str
        Name of object in field

    filter : str
        Filter name of images being stacked

    Note
    ----
    Produces a file with the following naming scheme ``object+'_'+filter+'_stack.list'``

    """

    # img, ota, filt = np.loadtxt('derived_props.txt', usecols=(0,1,2), dtype=str, unpack=True)
    # fwhm, zp_med, zp_std, bg_mean, bg_med, bg_std = np.loadtxt('derived_props.txt', usecols=(3,4,5,6,7,8), dtype=float, unpack=True)
    # 
    # keep = np.where(filt==filter)
    scaled_imgs = glob.glob(odi.scaledpath+'*'+filter+'*.fits')
    #sort glob list to match order in 'derived_props.txt'
    #this will sort on the dither number in the config file.
    scaled_imgs = sorted(scaled_imgs, key=lambda x: x[17:18].strip('_'))
    # head = scaled_imgs[0][:14]
    # Need to make a list of tails. The Job IDs will not always be the same
    # if different QR jobs were run (e.g. mix of user and operator images).
    # old definition tail = scaled_imgs[0][25:]
    # tail = []
    # for name in scaled_imgs:
    #     tail.append(name[25:])
    if not os.path.isfile(object.replace(' ','_')+'_'+filter+'_stack.list'):
        with open(object.replace(' ','_')+'_'+filter+'_stack.list','w+') as stack_file:
            for j,im in enumerate(scaled_imgs):
                # guide = odi.is_guide_ota(im, ota[j])
                # # print head+ota[j]+'.'+im+tail, factor
                # if not guide:
                print >> stack_file, im


def stack_images(stackname, refimg):
    """
    Stack the images that are in the list produced by ``make_stack_list`` using
    the ``IRAF`` task ``imcombine``. The following are the parameters used by
    ``imcombine``.

    - combine='average'
    - reject='none'
    - offsets='wcs'
    - masktype='goodvalue'
    - maskval=0
    - blank=-999
    - scale='none'
    - zero='none'
    - lthresh=-900
    - hthresh=60000
    - logfile=ota+'_stack.log'

    Parameters
    ----------
    stackname : str
        Name given to final stacked images
    refimg: str
        Name of reference image used in background calculation
    """
    from astropy.io import fits
    from pyraf import iraf
    print refimg
    fitsref = fits.open(refimg.f)
    hduref = fitsref[0]
    objname = stackname.replace(' ','_') #hduref.header['object'].replace(' ','_')
    filter_name = hduref.header['filter']
    ref_airmass = hduref.header['airmass']
    sky_med, sky_mean, sky_std = odi.find_new_bg(refimg, filter_name)
    odi.make_stack_list(objname, filter_name, refimg.inst)
    # sky_med = hduref.header['skybg']
    output = objname+'_'+filter_name+'.fits'
    output_bpm = objname+'_'+filter_name+'_bpm.pl'
    print '@'+objname+'_'+filter_name+'_stack.list'
    if not os.path.isfile(output):
        iraf.unlearn(iraf.immatch.imcombine, iraf.imutil.imarith)
        iraf.immatch.imcombine('@'+objname+'_'+filter_name+'_stack.list', 'temp', combine='average', reject='none', offsets='wcs', masktype='goodvalue', maskval=0, blank=-999, scale='none', zero='none', lthresh=-900, hthresh=60000)
        # iraf.imutil.imarith.setParam('operand1','temp')
        # iraf.imutil.imarith.setParam('op','+')
        # iraf.imutil.imarith.setParam('operand2',sky_med)
        # iraf.imutil.imarith.setParam('result',output)
        # iraf.imutil.imarith.setParam('verbose','yes')
        # iraf.imutil.imarith(mode='h')

        # flip the image so it's N-up E-left
        # first get the image dimensions from the header
        fitsstack = fits.open('temp.fits')
        xdim = fitsstack[0].header['NAXIS1']
        ydim = fitsstack[0].header['NAXIS2']
        iraf.imcopy('temp.fits['+repr(xdim)+':1,1:'+repr(ydim)+']', 'temp_flip.fits')

        iraf.imutil.imexpr('(a != -999) ? a + b : -999',output,'temp_flip.fits',sky_med)
        iraf.imutil.imexpr('a < 0',output_bpm, output)
        iraf.imutil.imdelete('temp, temp_flip', verify='no')
        iraf.unlearn(iraf.imutil.hedit)
        iraf.imutil.hedit.setParam('images',output)
        iraf.imutil.hedit.setParam('fields','BPM')
        iraf.imutil.hedit.setParam('value',output_bpm)
        iraf.imutil.hedit.setParam('add','yes')
        iraf.imutil.hedit.setParam('addonly','no')
        iraf.imutil.hedit.setParam('verify','no')
        iraf.imutil.hedit.setParam('update','yes')
        iraf.imutil.hedit(show='no', mode='h')

        # change the final CTYPENs to be TANs if they aren't already
        iraf.imutil.hedit.setParam('images',output)
        iraf.imutil.hedit.setParam('fields','CTYPE1')
        iraf.imutil.hedit.setParam('value','RA---TAN')
        iraf.imutil.hedit.setParam('add','yes')
        iraf.imutil.hedit.setParam('addonly','no')
        iraf.imutil.hedit.setParam('verify','no')
        iraf.imutil.hedit.setParam('update','yes')
        iraf.imutil.hedit(show='no', mode='h')

        iraf.imutil.hedit.setParam('images',output)
        iraf.imutil.hedit.setParam('fields','CTYPE2')
        iraf.imutil.hedit.setParam('value','DEC--TAN')
        iraf.imutil.hedit.setParam('add','yes')
        iraf.imutil.hedit.setParam('addonly','no')
        iraf.imutil.hedit.setParam('verify','no')
        iraf.imutil.hedit.setParam('update','yes')
        iraf.imutil.hedit(show='no', mode='h')

        # delete the inherited PV keywords
        # leaving them in will give you trouble with the stacked img wcs
        iraf.unlearn(iraf.imutil.hedit)
        iraf.imutil.hedit.setParam('images',output)
        iraf.imutil.hedit.setParam('fields','PV*')
        iraf.imutil.hedit.setParam('delete','yes')
        iraf.imutil.hedit.setParam('verify','no')
        iraf.imutil.hedit.setParam('update','yes')
        iraf.imutil.hedit(show='no', mode='h')

        # update the sky value and airmass keywords to match what they should be from the reference image
        iraf.unlearn(iraf.imutil.hedit)
        iraf.imutil.hedit.setParam('images',output)
        iraf.imutil.hedit.setParam('fields','airmass')
        iraf.imutil.hedit.setParam('value',ref_airmass)
        iraf.imutil.hedit.setParam('add','yes')
        iraf.imutil.hedit.setParam('addonly','no')
        iraf.imutil.hedit.setParam('verify','no')
        iraf.imutil.hedit.setParam('update','yes')
        iraf.imutil.hedit(show='no', mode='h')

        iraf.unlearn(iraf.imutil.hedit)
        iraf.imutil.hedit.setParam('images',output)
        iraf.imutil.hedit.setParam('fields','sky_medi')
        iraf.imutil.hedit.setParam('value',sky_med)
        iraf.imutil.hedit.setParam('add','yes')
        iraf.imutil.hedit.setParam('addonly','no')
        iraf.imutil.hedit.setParam('verify','no')
        iraf.imutil.hedit.setParam('update','yes')
        iraf.imutil.hedit(show='no', mode='h')

        iraf.unlearn(iraf.imutil.hedit)
        iraf.imutil.hedit.setParam('images',output)
        iraf.imutil.hedit.setParam('fields','sky_mean')
        iraf.imutil.hedit.setParam('value',sky_mean)
        iraf.imutil.hedit.setParam('add','yes')
        iraf.imutil.hedit.setParam('addonly','no')
        iraf.imutil.hedit.setParam('verify','no')
        iraf.imutil.hedit.setParam('update','yes')
        iraf.imutil.hedit(show='no', mode='h')

        iraf.unlearn(iraf.imutil.hedit)
        iraf.imutil.hedit.setParam('images',output)
        iraf.imutil.hedit.setParam('fields','sky_std')
        iraf.imutil.hedit.setParam('value',sky_std)
        iraf.imutil.hedit.setParam('add','yes')
        iraf.imutil.hedit.setParam('addonly','no')
        iraf.imutil.hedit.setParam('verify','no')
        iraf.imutil.hedit.setParam('update','yes')
        iraf.imutil.hedit(show='no', mode='h')

        # finally reset the LTV, LTM keywords so that the physical coordinate mapping is correct
        iraf.unlearn(iraf.imutil.hedit)
        iraf.imutil.hedit.setParam('images',output)
        iraf.imutil.hedit.setParam('fields','LTV1')
        iraf.imutil.hedit.setParam('value',0.)
        iraf.imutil.hedit.setParam('add','yes')
        iraf.imutil.hedit.setParam('addonly','no')
        iraf.imutil.hedit.setParam('verify','no')
        iraf.imutil.hedit.setParam('update','yes')
        iraf.imutil.hedit(show='no', mode='h')

        iraf.unlearn(iraf.imutil.hedit)
        iraf.imutil.hedit.setParam('images',output)
        iraf.imutil.hedit.setParam('fields','LTM1_1')
        iraf.imutil.hedit.setParam('value',1.)
        iraf.imutil.hedit.setParam('add','yes')
        iraf.imutil.hedit.setParam('addonly','no')
        iraf.imutil.hedit.setParam('verify','no')
        iraf.imutil.hedit.setParam('update','yes')
        iraf.imutil.hedit(show='no', mode='h')

        iraf.unlearn(iraf.imutil.hedit)
        iraf.imutil.hedit.setParam('images',output)
        iraf.imutil.hedit.setParam('fields','LTM2_2')
        iraf.imutil.hedit.setParam('value',1.)
        iraf.imutil.hedit.setParam('add','yes')
        iraf.imutil.hedit.setParam('addonly','no')
        iraf.imutil.hedit.setParam('verify','no')
        iraf.imutil.hedit.setParam('update','yes')
        iraf.imutil.hedit(show='no', mode='h')

    return output

def instrument(img):
    """
    A function to grab what version of ODI has been used.

    Parameters
    ----------
    img : str
        Name of image

    Returns
    -------
    instrument_name : str
        Will return ``5odi`` or ``podi``.
    """

    from astropy.io import fits
    hdulist = fits.open(img.f)
    instrument_name = hdulist[0].header['INSTRUME']
    hdulist.close()
    print 'Setting instrument to: ', instrument_name

    if instrument_name == '5odi':
        # odi.OTA_dictionary = odi.odi5narrow_dictionary
        odi.OTA_dictionary = odi.odi5_dictionary
    else :
        odi.OTA_dictionary = odi.podi_dictionary

    return instrument_name

def qr_img_lists(odi_filter):
    """
    (Currently not used)
    a function to build the image lists that will be fed into
    odi_process.py and other routines. This will be based on
    the filters provided.

    Filter example = odi_g

    """
    glob_cmd = '*_'+odi_filter+'*.fits'
    img_list = glob.glob(glob_cmd)
    img_list.sort()

    return img_list

def tan_header_fix(hdu):
    if 'TAN' in hdu.header['CTYPE1']:
        pvlist = hdu.header['PV*']
        for pv in pvlist:
            hdu.header.delete_keyword(pv)
    return hdu

def tpv2tan_hdr(img, ota):
    image = odi.reprojpath+'reproj_'+ota+'.'+img.stem()
    # change the CTYPENs to be TANs if they aren't already
    print 'TPV -> TAN in ', image
    iraf.imutil.hedit.setParam('images',image)
    iraf.imutil.hedit.setParam('fields','CTYPE1')
    iraf.imutil.hedit.setParam('value','RA---TAN')
    iraf.imutil.hedit.setParam('add','yes')
    iraf.imutil.hedit.setParam('addonly','no')
    iraf.imutil.hedit.setParam('verify','no')
    iraf.imutil.hedit.setParam('update','yes')
    iraf.imutil.hedit(show='no', mode='h')

    iraf.imutil.hedit.setParam('images',image)
    iraf.imutil.hedit.setParam('fields','CTYPE2')
    iraf.imutil.hedit.setParam('value','DEC--TAN')
    iraf.imutil.hedit.setParam('add','yes')
    iraf.imutil.hedit.setParam('addonly','no')
    iraf.imutil.hedit.setParam('verify','no')
    iraf.imutil.hedit.setParam('update','yes')
    iraf.imutil.hedit(show='no', mode='h')

    # delete any PV keywords
    # leaving them in will give you trouble with the img wcs
    iraf.unlearn(iraf.imutil.hedit)
    iraf.imutil.hedit.setParam('images',image)
    iraf.imutil.hedit.setParam('fields','PV*')
    iraf.imutil.hedit.setParam('delete','yes')
    iraf.imutil.hedit.setParam('verify','no')
    iraf.imutil.hedit.setParam('update','yes')
    iraf.imutil.hedit(show='no', mode='h')

def find_ref_image(images):
    imgs, fwhm, zp_med, zp_std, bg_mean, bg_median, bg_std = np.loadtxt('derived_props.txt', usecols=(0,3,4,5,6,7,8), unpack=True)
    filter_string = np.loadtxt('derived_props.txt', usecols=(2,), unpack=True,dtype=str)

    lvls = []
    ams = []
    zps = []
    #print images
    print '#'+ ' '*len(images[0].f)+'     bg         airmass'
    for j,im in enumerate(images):
        hdulist = odi.fits.open(im.f)
        airmass = hdulist[0].header['AIRMASS']
        filter  = hdulist[0].header['FILTER']
        these = np.where((imgs.astype(int)==int(im.dither())) & (filter_string == filter))
        bg_lvl = np.mean(bg_median[these])
        lvls.append(bg_lvl)
        ams.append(airmass)
        hdulist.close()
        print im.dither(), im.f, '%10.3f'%bg_lvl, '%10.3f'%airmass
        ref_img = np.argmin(np.array(ams))
    print 'reference image:',images[ref_img].stem()
    return ref_img
    
def imalign(images, square=False):
    from astropy.wcs import WCS
    from astropy.io import fits

    img_ref = images[0]
    img     = images[1]

    # fetch a gaia catalog for the full image footprint
    outputg = img_ref.f+'.gaia'
    gaia_cat = odi.get_gaia_coords(images[0], ota='None', inst='podi', output=outputg)
    # print gaia_cat
    # convert the ra, dec to image coordinates in each image
    # first get the wcs for each image
    # (pick the first image as a "reference")
    hdu_ref = fits.open(img_ref.f)
    naxis1_ref = hdu_ref[0].header['NAXIS1']
    naxis2_ref = hdu_ref[0].header['NAXIS2']
    w_ref = WCS(hdu_ref[0].header)
    x_ref, y_ref = w_ref.all_world2pix(gaia_cat.ra, gaia_cat.dec, 1)

    hdu_img = fits.open(img.f)
    naxis1 = hdu_img[0].header['NAXIS1']
    naxis2 = hdu_img[0].header['NAXIS2']
    w_img = WCS(hdu_img[0].header)
    x_img, y_img = w_img.all_world2pix(gaia_cat.ra, gaia_cat.dec, 1)

    # compute the pair-wise integer pixel shifts between the image and the reference
    x_shift, y_shift = np.rint(np.median(x_ref-x_img)), np.rint(np.median(y_ref-y_img))

    # shift the images so that they are aligned to the "reference"
    # use relative coordinates-- negative values are applied to the image, positives to the 'reference'
    iraf.imcopy(img.f, 'temp.fits')
    iraf.imcopy(img_ref.f, 'temp_ref.fits')

    if x_shift < 0.0:
        trim_img = 'temp.fits[{:d}:{:d},*]'.format(int(abs(x_shift))+1, naxis1)
        iraf.imcopy(trim_img, 'temp.fits')
    else:
        trim_img = 'temp_ref.fits[{:d}:{:d},*]'.format(int(abs(x_shift))+1, naxis1_ref)
        iraf.imcopy(trim_img, 'temp_ref.fits')

    if y_shift < 0.0:
        trim_img = 'temp.fits[*,{:d}:{:d}]'.format(int(abs(y_shift))+1, naxis2)
        iraf.imcopy(trim_img, 'temp.fits')
    else:
        trim_img = 'temp_ref.fits[*,{:d}:{:d}]'.format(int(abs(y_shift))+1, naxis2_ref)
        iraf.imcopy(trim_img, 'temp_ref.fits')

    # then take any excess pixels off at the high end of the range in each dimension so that the images have identical dimensions
    # select the smallest dimension in any image
    new_img = img.f[:-5]+'_match.fits'
    new_ref = img_ref.f[:-5]+'_match.fits'

    hdu_tempref = fits.open('temp_ref.fits')
    naxis1_tempref = hdu_tempref[0].header['NAXIS1']
    naxis2_tempref = hdu_tempref[0].header['NAXIS2']

    hdu_temp = fits.open('temp.fits')
    naxis1_temp = hdu_temp[0].header['NAXIS1']
    naxis2_temp = hdu_temp[0].header['NAXIS2']

    min_xsize = min(naxis1_temp, naxis1_tempref)
    min_ysize = min(naxis2_temp, naxis2_tempref)

    trim_img = 'temp.fits[1:{:d},1:{:d}]'.format(int(min_xsize)-1, int(min_ysize)-1) 
    trim_ref = 'temp_ref.fits[1:{:d},1:{:d}]'.format(int(min_xsize)-1, int(min_ysize)-1)

    iraf.imcopy(trim_img, new_img)
    iraf.imcopy(trim_ref, new_ref)    

    # delete the temporary working images
    iraf.imdelete('temp.fits')
    iraf.imdelete('temp_ref.fits')
    
def main():
    images = [odi.StackedImage('AGC198511_odi_g.fits'),odi.StackedImage('AGC198511_odi_i.fits')]
    imalign(images)

if __name__ == '__main__':
    main()
