#!/usr/bin/env python
import sys, os, glob, string
import numpy as np
import astropy as ast
import matplotlib.pyplot as plt
from pyraf import iraf
from tqdm import tqdm
import odi_config as odi

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
    hdulist = odi.fits.open(img.f)
    existing_radesys = hdulist[0].header['RADESYS']
    print img.f
    print '--> Existing RADESYS value:', existing_radesys
    for k in img.otas.keys():
        existing_ctype1 = hdulist[k].header['CTYPE1']
        existing_ctype2 = hdulist[k].header['CTYPE2']
        print '--> Existing CTYPE1 value:', existing_ctype1
        print '--> Existing CTYPE2 value:', existing_ctype2

def main():
    object_str, filters, instrument, images, illcor_flag, skyflat_src, wcs_flag, reproject_flag, scale_flag, stack_flag, gaia_flag, cluster_flag, ra_center, dec_center, min_radius = odi.cfgparse('config.yaml', verbose=False)
    
    images_ = [img for sublist in images.values() for img in sublist]
    for img  in images_:
        repair_wcs_keywords(img)

if __name__ == '__main__':
    main()
