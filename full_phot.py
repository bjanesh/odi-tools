import sys, os, glob, string
import numpy as np
import astropy as ast
import matplotlib.pyplot as plt
from pyraf import iraf
import odi_config as odi
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
from collections import OrderedDict

def find_sources_full(img,fwhm,bg_std,threshold=4.0):
    """
    Use ``pyraf daofind`` to located sources on a stacked image.
    ``doafind`` options ::
    iraf.unlearn(iraf.apphot.daofind)
    iraf.datapars.setParam('fwhmpsf',fwhm,check=1)
    iraf.datapars.setParam('datamin',-900,check=1)
    iraf.datapars.setParam('datamax',60000,check=1)
    iraf.datapars.setParam('sigma',bg_std,check=1)
    iraf.findpars.setParam('threshold',threshold)
    iraf.apphot.daofind.setParam('output',output)
    iraf.apphot.daofind(image=img, verbose="no", verify='no')

    Parameters
    ----------
    img : str
         String containing name of the image currently in use

    fwhm : float
        fwhm measure of sources in field

    bg_std : float
        standard deviation of background in image

    threshold : float
        detection threshold for sources

    Note
    ----
    Produces a coordinate file based on the name of the image.
    The file name will be ``img.nofits()+_sources.coo``
    """
    output = img.nofits()+'_sources.coo'
    if not os.path.isfile(output):
        print('Locating sources on ',img)
        print('Will output to ',output)
        iraf.unlearn(iraf.apphot.daofind)
        iraf.datapars.setParam('fwhmpsf',fwhm,check=1)
        iraf.datapars.setParam('datamin',-900,check=1)
        iraf.datapars.setParam('datamax',60000,check=1)
        iraf.datapars.setParam('sigma',bg_std,check=1)
        iraf.findpars.setParam('threshold',threshold)
        iraf.apphot.daofind.setParam('output',output)
        iraf.apphot.daofind(image=img, verbose="no", verify='no')

def phot_sources_full(img,fwhm,airmass,apfactor):
    """
    Run ``pyraf phot`` on the sources found by ``find sources_full``

    Parameters
    ----------

    img : str
        String containing name of the image currently in use

    fwhm : float
        fwhm measurement in image

    airmass : float
        airmass of image

    apfactor : float
        multile of fwhm to use for photometry

    Note
    ----
    Will retrun a ``.phot`` table with the name ``img.nofits()+.srcphot``

    """
    iraf.ptools(_doprint=0)
    coords = img.nofits()+'_sources.coo'
    output = img.nofits()+'.phot.1'
    phot_tbl = img.nofits()+'.srcphot'
    if not os.path.isfile(phot_tbl) :
        print('phot-ing ', img, ' from daofind')
        iraf.unlearn(iraf.apphot.phot,iraf.datapars,iraf.photpars,iraf.centerpars,iraf.fitskypars)
        iraf.apphot.phot.setParam('interactive',"no")
        iraf.apphot.phot.setParam('verify',"no")
        iraf.datapars.setParam('datamax',50000.)
        iraf.datapars.setParam('gain',"gain")
        iraf.datapars.setParam('ccdread','rdnoise')
        iraf.datapars.setParam('exposure',"exptime")
        iraf.datapars.setParam('xairmass',airmass)
        iraf.datapars.setParam('filter',"filter")
        iraf.datapars.setParam('obstime',"time-obs")
        iraf.datapars.setParam('sigma',"INDEF")
        iraf.photpars.setParam('zmag',0.)
        iraf.centerpars.setParam('cbox',9.)
        iraf.centerpars.setParam('maxshift',3.)
        iraf.fitskypars.setParam('salgorithm',"median")
        iraf.fitskypars.setParam('dannulus',10.)

        iraf.datapars.setParam('fwhmpsf',fwhm)
        iraf.photpars.setParam('apertures',apfactor*fwhm)
        iraf.fitskypars.setParam('annulus',6.5*fwhm)

        iraf.apphot.phot(image=img, coords=coords, output=output)

        with open(phot_tbl,'w+') as txdump_out:
            iraf.ptools.txdump(textfiles=output, fields="id,mag,merr,msky,stdev,rapert,xcen,ycen,ifilter,xairmass,image",expr='yes', headers='no', Stdout=txdump_out)

        outputfile_clean = open(phot_tbl.replace('.srcphot','_clean.srcphot'),"w")
        for line in open(phot_tbl,"r"):
            if not 'INDEF' in line:
                outputfile_clean.write(line)
            if 'INDEF' in line:
                outputfile_clean.write(line.replace('INDEF','999'))
        outputfile_clean.close()
        os.rename(phot_tbl.replace('.srcphot','_clean.srcphot'),phot_tbl)

def phot_sources_xy2sky(img,inst):
    """
    Convert the x,y positions in the phot table produced by
    ``phot_sources_full`` into Ra and Dec positions.

    Parameters
    ----------

    img : str
        String containing name of the image currently in use

    inst : str
        ODI configuration. ``podi`` or ``5odi``

    Note
    ----
    Returns a table with the name ``img.nofits()+.srcphotrd``

    """
    phot_tbl = img.nofits()+'.srcphot'
    outputradec = img.nofits()+'.srcphotrd'
    hdulist= odi.fits.open(img.f)

    #if inst == 'podi':
    #    header = hdulist[0].header
    #    print header['CTYPE1']
    #    print header['CTYPE2']
    #    header['CTYPE1'] = 'RA---TPV'
    #    header['CTYPE2'] = 'DEC--TPV'
    #    print header['CTYPE1']
    #    print header['CTYPE2']
    #    w = odi.WCS(header)
        #w.wcs.ctype = ["RA---TPV", "DEC--TPV"]
    #     print w
        #pvlist = hdulist[0].header['PV*']
        #for pv in pvlist:
            #tpv = 'T'+pv
            #hdulist[0].header.rename_keyword(pv, tpv, force=False)
    w = odi.WCS(hdulist[0].header)
    MAG, MERR, SKY, SERR, RAPERT, XPOS, YPOS = np.loadtxt(phot_tbl, usecols=(1,2,3,4,5,6,7), dtype=float, unpack=True)
    with open(outputradec, 'w+') as fxy:
        for i,c in enumerate(XPOS):
            coords2 = [[XPOS[i],YPOS[i]]]
            pixcrd2 = w.wcs_pix2world(coords2, 1)
            print(pixcrd2[0][0], pixcrd2[0][1],XPOS[i],YPOS[i],MAG[i], MERR[i],SKY[i],SERR[i],RAPERT[i], file=fxy)
    hdulist.close()

def match_phot_srcs(img1,img2):
    """
    Match the sources in two images. This function reads in the photometry tables
    produces by ``phot_sources_xy2sky`` and used the Ra and Dec positions to
    match the sources between the images.

    Parameters
    ----------
    img1 : str
        Name of the stacked image in the first filter (e.g. odi_g)
    img2 : str
        Name of the stacked image in the second filter (e.g. odi_r)

    Note
    ----
    Produces a catalog of matched sources for each image.
    ``img1[:-5]+.match.srscrd`` and ``img2[:-5]+.match.srscrd`` These magnitudes
    are combined in the file ``calibration.dat``.

    """

    img1_srcs =img1[:-5]+'.srcphotrd'
    img1_srsc_match = img1[:-5]+'.match.srscrd'
    img2_srcs =img2[:-5]+'.srcphotrd'
    img2_srsc_match = img2[:-5]+'.match.srscrd'

    ra_1, dec_1,x_1,y_1,mag_1,merr_1,sky_1,serr_1,rapert_1 = np.loadtxt(img1_srcs,usecols=(0,1,2,3,4,5,6,7,8),unpack=True)
    ra_2, dec_2,x_2,y_2,mag_2,merr_2,sky_2,serr_2,rapert_2 = np.loadtxt(img2_srcs,usecols=(0,1,2,3,4,5,6,7,8),unpack=True)


    img1_catalog = SkyCoord(ra = ra_1*u.degree, dec= dec_1*u.degree)

    img2_catalog = SkyCoord(ra = ra_2*u.degree, dec= dec_2*u.degree)

    id_img1, id_img2, d2d, d3d = img2_catalog.search_around_sky(img1_catalog,0.00005*u.deg)

    print(len(id_img1),len(id_img2),len(x_1),len(x_2))

    ra_1      =    ra_1[id_img1]
    dec_1     =    dec_1[id_img1]
    x_1       =    x_1[id_img1]
    y_1       =    y_1[id_img1]
    mag_1     =    mag_1[id_img1]
    merr_1    =    merr_1[id_img1]
    sky_1     =    sky_1[id_img1]
    serr_1    =    serr_1[id_img1]
    rapert_1  =    rapert_1[id_img1]

    ra_2      =    ra_2[id_img2]
    dec_2     =    dec_2[id_img2]
    x_2       =    x_2[id_img2]
    y_2       =    y_2[id_img2]
    mag_2     =    mag_2[id_img2]
    merr_2    =    merr_2[id_img2]
    sky_2     =    sky_2[id_img2]
    serr_2    =    serr_2[id_img2]
    rapert_2  =    rapert_2[id_img2]

    with open(img1_srsc_match,'w+') as m1:
        with open(img2_srsc_match,'w+') as m2:
            with open('calibration.dat','w+') as cal:
                for i,s in enumerate(ra_1):
                    print(ra_1[i], dec_1[i],x_1[i],y_1[i],mag_1[i],merr_1[i],sky_1[i],serr_1[i],rapert_1[i], file=m1)
                    print(ra_2[i], dec_2[i],x_2[i],y_2[i],mag_2[i],merr_2[i],sky_2[i],serr_2[i],rapert_2[i], file=m2)
                    print(i, x_1[i], y_1[i], rapert_1[i], mag_1[i], merr_1[i], x_2[i], y_2[i], rapert_2[i], mag_2[i], merr_2[i], file=cal)
    m1.close()
    m2.close()
    junk = open('match.reg',"w")
    for i in range(len(ra_1)):
        print(x_2[i], y_2[i], file=junk)
    junk.close()

def calc_calibrated_mags(apcor_g, cal_A_g, apcor_r, cal_A_r, photcalFile, object_str):
    """
    Convert the instrumental magnitudes to calibrated magnitudes using the help
    file produce by ``full_calibrate``.

    Parameters
    ----------

    apcor_g : float
        aperture correction in the first filter

    cal_A_g : float
        extinction in first filter

    apcor_r : float
        aperture correction in the second filter

    cal_A_r : float
        extinction in second filter

    photcalFile : str
        help file produced by ``full_calibrate``

    object_str : str
        object string (e.g. GCPair-F1)


    Note
    ----
    The calibrated magnitudes in both filters will be save to a file called
    ``calibrated_mags.dat``.


    """
    from matplotlib import gridspec
    kg = 0.200
    kr = 0.12
    ki = 0.058
    rdnoise = 6.5

    # get the photometric calibration coefficients from Steven's help file <--
    # or from the image header/fits table/ whatever
    photcalFile = open(photcalFile)
    photcal = photcalFile.read()
    photcalLines = photcal.splitlines()

    mu_gr = float(photcalLines[28].split()[5])
    zp_gr = float(photcalLines[30].split()[4])
    eps_gr = float(photcalLines[34].split()[5])
    zp_r = float(photcalLines[36].split()[4])
    amg = float(photcalLines[25].split()[5])
    amr = float(photcalLines[26].split()[5])
    photcalFile.close()

    nid,gx,gy,g_i,g_ierr,rx,ry,r_i,r_ierr = np.loadtxt('calibration.dat',usecols=(0,1,2,4,5,6,7,9,10),unpack=True)
    print(apcor_g,apcor_r)
    g0 = g_i - (kg*amg) + apcor_g
    r0 = r_i - (kr*amr) + apcor_r
    gmr = mu_gr*(g0-r0) + zp_gr

    r_mag = r0 + eps_gr*gmr + zp_r
    g_mag = gmr + r_mag - cal_A_g
    r_mag = r_mag - cal_A_r
    gmr = g_mag - r_mag

    print('Median (g-r) :: g - r = {0:7.4f}'.format(np.median(gmr)))
    print('Final number of phot-ed stars :: g = {0:5d} : r = {1:5d}'.format(len(g_mag),len(r_mag)))

    g_mag_lims = [g_mag[i] for i in range(len(g_mag)) if (g_ierr[i] >= 0.2)]
    r_mag_lims = [r_mag[i] for i in range(len(r_mag)) if (r_ierr[i] >= 0.2)]
    with open('calibrated_magstest.dat', 'w+') as f3:
        for i in range(len(rx)) :
            print('{0:8.2f} {1:8.2f} {2:12.3f} {3:12.3f} {4:8.2f} {5:8.2f} {6:12.3f} {7:12.3f} {8:12.3f} '.format(gx[i],gy[i],g_mag[i],g_ierr[i],rx[i],ry[i],r_mag[i],r_ierr[i],gmr[i]), file=f3)

    plt.clf()

    fig = plt.figure(figsize=(10, 8))
    fig.subplots_adjust(hspace=0)
    gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])
    ax0 = plt.subplot(gs[0])

    ax0.scatter(gmr, r_mag, s=2, color='black', marker='o', edgecolors='none')
    ax0.set_ylabel('$r$')
    ax0.set_xlabel('$(g-r)$')
    ax0.set_ylim(24,10)
    ax0.set_xlim(-1,2)

    ax1 = plt.subplot(gs[1])

    ax1.scatter(r_ierr, r_mag, s=2, color='black', marker='o', edgecolors='none')
    # ax1.set_ylabel('$r$')
    ax1.set_xlabel('inst $r$ err.')
    ax1.set_ylim(24,10)
    ax1.set_xlim(-0.002,0.05)
    plt.setp(ax1.get_yticklabels(), visible=False)

    plt.tight_layout()
    plt.savefig(object_str+'_CMD.pdf')
