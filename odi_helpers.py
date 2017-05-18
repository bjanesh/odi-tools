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

    img, ota, filt = np.loadtxt('derived_props.txt', usecols=(0,1,2), dtype=str, unpack=True)
    fwhm, zp_med, zp_std, bg_mean, bg_med, bg_std = np.loadtxt('derived_props.txt', usecols=(3,4,5,6,7,8), dtype=float, unpack=True)

    keep = np.where(filt==filter)
    bg_std_med, bg_std_std = np.median(bg_std[keep]), np.std(bg_std[keep])
    # print bg_std_med, bg_std_std
    scaled_imgs = glob.glob(odi.scaledpath+'*'+filter+'*.fits')
    #sort glob list to match order in 'derived_props.txt'
    #this will sort on the dither number in the config file.
    scaled_imgs = sorted(scaled_imgs, key=lambda x: int(x[24]))
    head = scaled_imgs[0][:14]
    # Need to make a list of tails. The Job IDs will not always be the same
    # if different QR jobs were run (e.g. mix of user and operator images).
    # old definition tail = scaled_imgs[0][25:]
    tail = []
    for name in scaled_imgs:
        tail.append(name[25:])
    if not os.path.isfile(object+'_'+filter+'_stack.list'):
        with open(object+'_'+filter+'_stack.list','w+') as stack_file:
            if inst == 'podi':
                # don't exclude any otas from podi data, the guide otas are probably not here anyway
                for j,im in enumerate(img[keep]):
                    print >> stack_file, head+ota[j]+'.'+im+tail[j]
            else:
                for j,im in enumerate(img[keep]):
                    factor = np.absolute(bg_std[keep][j] - bg_std_med)/bg_std_std
                    # print head+ota[j]+'.'+im+tail, factor
                    if factor < 2.:
                        print >> stack_file, head+ota[j]+'.'+im+tail[j]


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
    objname = hduref.header['object'].replace(' ','_')
    filter_name = hduref.header['filter']
    ref_airmass = hduref.header['airmass']
    sky_med, sky_mean, sky_std = odi.find_new_bg(refimg, filter_name)
    odi.make_stack_list(objname, filter_name, refimg.inst)
    # sky_med = hduref.header['skybg']
    output = stackname+'_'+filter_name+'.fits'
    output_bpm = stackname+'_'+filter_name+'_bpm.pl'
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
    # fetch a gaia catalog for the full image footprint
    outputg = images[0].f+'.gaia'
    gaia_cat = odi.get_gaia_coords(images[0], ota='None', inst='5odi', output=outputg)
    # print gaia_cat
    # convert the ra, dec to image coordinates in each image
    # first get the wcs for each image
    # (pick the first image as a "reference")
    hdu_ref = fits.open(images[0].f)
    w_ref = WCS(hdu_ref[0].header)
    x_ref, y_ref = w_ref.all_world2pix(gaia_cat.ra, gaia_cat.dec, 1)
    # compute the pair-wise integer pixel shifts for each image
    x_shift = {}
    y_shift = {}
    naxis1 = {}
    naxis2 = {}
    for img in images:
        hdu_img = fits.open(img.f)
        naxis1[img.f] = hdu_img[0].header['NAXIS1']
        naxis2[img.f] = hdu_img[0].header['NAXIS2']
        w_img = WCS(hdu_img[0].header)
        x_img, y_img = w_img.all_world2pix(gaia_cat.ra, gaia_cat.dec, 1)
        x_shift[img.f], y_shift[img.f] = np.rint(np.median(x_ref-x_img)), np.rint(np.median(y_ref-y_img))
    
    
    print x_shift, y_shift, naxis1, naxis2
    # shift the images so that they are aligned to the "reference"
    max_xshift = max(np.abs(x_shift.values()))
    max_yshift = max(np.abs(y_shift.values()))
    new_x_hdim = int(naxis1[images[0].f] - 2*max_xshift - 2)/2
    new_y_hdim = int(naxis2[images[0].f] - 2*max_yshift - 2)/2
    
    for img in images:
        xcent, ycent = int((naxis1[img.f]/2)-x_shift[img.f]), int((naxis2[img.f]/2)-y_shift[img.f])
        if square: 
            ledge, redge = xcent-min(new_x_hdim, new_y_hdim), xcent+min(new_x_hdim, new_y_hdim)
            bedge, tedge = ycent-min(new_x_hdim, new_y_hdim), ycent+min(new_x_hdim, new_y_hdim)
        else:
            ledge, redge = xcent-new_x_hdim, xcent+new_x_hdim
            bedge, tedge = ycent-new_y_hdim, ycent+new_y_hdim
        trim_img = '{:s}[{:d}:{:d},{:d}:{:d}]'.format(img.f, ledge, redge, bedge, tedge)
        new_img = img.f[:-5]+'_match.fits'
        iraf.imcopy(trim_img, new_img)
    
def main():
    images = [odi.StackedImage('AGC198511_odi_g.fits'),odi.StackedImage('AGC198511_odi_i.fits')]
    imalign(images)

if __name__ == '__main__':
    main()
