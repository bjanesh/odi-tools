import sys, os, glob, string
import numpy as np
import astropy as ast
import matplotlib.pyplot as plt
from pyraf import iraf
from tqdm import tqdm
import odi_config as odi

def zeropoint_ota(img, ota, fwhm):
    iraf.ptools(_doprint=0)
    # otaext = {'33':1,'34':2,'44':3,'43':4,'42':5,'32':6,'22':7,'23':8,'24':9}
    # values determined by ralf/daniel @ wiyn
    kg = 0.20
    kr = 0.12
    ki = 0.058

    # first grab the header and hang on to it so we can use other values
    # hdulist = ast.io.fits.open(img)
    # hdr = hdulist[ota].header
    # hdulist.close()

    image = odi.reprojpath+'reproj_'+ota+'.'+img.stem()
    coords = odi.coordspath+'reproj_'+ota+'.'+img.base()+'.sdssxy'
    output = odi.coordspath+img.nofits()+'.'+ota+'.phot.1'
    phot_tbl = odi.coordspath+img.nofits()+'.'+ota+'.sdssphot'

    # alas, we must use IRAF apphot to do the measuring
    # first set common parameters (these shouldn't change if you're using ODI)
    iraf.unlearn(iraf.phot,iraf.datapars,iraf.photpars,iraf.centerpars,iraf.fitskypars)
    iraf.apphot.phot.setParam('interactive',"no")
    iraf.apphot.phot.setParam('verify',"no")
    iraf.datapars.setParam('datamax',50000.)
    iraf.datapars.setParam('gain',"gain")
    iraf.datapars.setParam('ccdread',"rdnoise")
    iraf.datapars.setParam('exposure',"exptime")

    iraf.datapars.setParam('filter',"filter")
    iraf.datapars.setParam('obstime',"time-obs")
    iraf.datapars.setParam('sigma',"INDEF")
    iraf.photpars.setParam('zmag',0.)
    iraf.centerpars.setParam('cbox',9.)
    iraf.centerpars.setParam('maxshift',3.)
    iraf.fitskypars.setParam('salgorithm',"median")
    iraf.fitskypars.setParam('dannulus',10.)

    iraf.datapars.setParam('airmass','airmass')
    iraf.datapars.setParam('fwhmpsf',fwhm)
    iraf.photpars.setParam('apertures',3.*fwhm) # use a big aperture for this
    iraf.fitskypars.setParam('annulus',4.*fwhm)

    if not os.path.isfile(output):
        iraf.apphot.phot(image=image, coords=coords, output=output)
    with open(phot_tbl,'w+') as txdump_out :
        iraf.ptools.txdump(textfiles=output, fields="id,mag,merr,msky,stdev,rapert,xcen,ycen,ifilter,xairmass,peak,flux,image", expr='yes', headers='no', Stdout=txdump_out)
    
    outputfile_clean = open(phot_tbl.replace('.sdssphot','_clean.sdssphot'),"w")
    for line in open(phot_tbl,"r"):
        if not 'INDEF' in line:
            outputfile_clean.write(line)
        if 'INDEF' in line:
            outputfile_clean.write(line.replace('INDEF','999'))
    outputfile_clean.close()
    os.rename(phot_tbl.replace('.sdssphot','_clean.sdssphot'),phot_tbl)
    
    gMAG, gMERR, gSKY, gSERR, gRAPERT, gXPOS, gYPOS = np.loadtxt(phot_tbl, usecols=(1,2,3,4,5,6,7), dtype=float, unpack=True)
    gXAIRMASS = np.loadtxt(phot_tbl, usecols=(9,), dtype=float, unpack=True)
    gFILTER = np.loadtxt(phot_tbl, usecols=(8,), dtype=str, unpack=True)
    gID = np.loadtxt(phot_tbl, usecols=(0,), dtype=int, unpack=True)

    # keep the actual ID number to select from SDSS stars
    keep = gID - 1

    # read in the catalog values
    x, y, ra, dec, u, ue, g, ge, r, re, i, ie, z, ze = np.loadtxt(coords, usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13), unpack=True)

    # pick out the ones that match the good phot stars
    g, ge, r, re, i, ie = np.array(g[keep]), np.array(ge[keep]), np.array(r[keep]), np.array(re[keep]), np.array(i[keep]), np.array(ie[keep])

    # and reduce the other vectors
    # gXPOS, gYPOS, gMAG, gMERR, gSKY, gSERR = np.array(gXPOS[keep]), np.array(gYPOS[keep]), np.array(gMAG[keep]), np.array(gMERR[keep]), np.array(gSKY[keep]), np.array(gSERR[keep])

    gXAIRMASS = gXAIRMASS[0]
    gRAPERT = gRAPERT[0]

    #podicut, sdsscut = 0.01, 0.03
    podicut, sdsscut = 0.01, 0.10
    # pick the right airmass extinction coefficient
    kdict = {'odi_g':kg, 'odi_r':kr, 'odi_i':ki}
    mdict = {'odi_g':g, 'odi_r':r, 'odi_i':i}
    k = kdict[gFILTER[0]]
    mag_cat = mdict[gFILTER[0]]
    g0 = gMAG - k*gXAIRMASS
    zp = mag_cat - g0
    color = g - r
    errcut = [j for j in range(len(gMERR)) if (gMERR[j] < podicut and ge[j] < sdsscut)]
    if len(g0[errcut]) < 5:
        print 'Too few stars for fit..'
        return 999.999, 999.999, phot_tbl
    else:
        print 'fitting wtih '+repr(len(g0[errcut]))+' stars...'
        
        # fit zero point
        # linear lsq with numpy.polyfit
        p, pcov = np.polyfit(color[errcut], zp[errcut], 1, cov=True)
        perr = np.sqrt(np.diag(pcov))
        eps_g, zp_g, std_eps_g, std_zp_g = p[0], p[1], perr[0], perr[1]
        
        print '--------------------------------------------------------------------------'
        print 'Here are the fit values:'
        print 'eps  '+'      std_eps  '+'  zp  '+'         std_zp'
        print '{0:10.7f} {1:10.7f} {2:10.7f} {3:10.7f}'.format(eps_g, std_eps_g, zp_g, std_zp_g)
        star_zp_g = zp - eps_g*color
        print 'std. dev. in ZP per star (not fit): {0:10.7f}'.format(np.std(star_zp_g[errcut]))
        # for i in range(len(zp)):
        #     print color[i], zp[i], gMERR[i], ge[i]
        
        print np.median(star_zp_g[errcut]), np.std(star_zp_g[errcut])
        return np.median(star_zp_g[errcut]), np.std(star_zp_g[errcut]), phot_tbl

def zeropoint_full(img, fwhm):
    iraf.ptools(_doprint=0)
    # otaext = {'33':1,'34':2,'44':3,'43':4,'42':5,'32':6,'22':7,'23':8,'24':9}
    # values determined by ralf/daniel @ wiyn
    kg = 0.20
    kr = 0.12
    ki = 0.058
    
    # first grab the header and hang on to it so we can use other values
    # hdulist = ast.io.fits.open(img)
    # hdr = hdulist[ota].header
    # hdulist.close()
    
    # image = odi.reprojpath+'reproj_'+ota+'.'+img.stem()
    coords = img.nofits()+'.sdssxy'
    output = img.nofits()+'.phot.1'
    phot_tbl = img.nofits()+'.sdssphot'
    
    # alas, we must use IRAF apphot to do the measuring
    # first set common parameters (these shouldn't change if you're using ODI)
    iraf.unlearn(iraf.phot,iraf.datapars,iraf.photpars,iraf.centerpars,iraf.fitskypars)
    iraf.apphot.phot.setParam('interactive',"no")
    iraf.apphot.phot.setParam('verify',"no")
    iraf.datapars.setParam('datamax',50000.)
    iraf.datapars.setParam('gain',"gain")
    iraf.datapars.setParam('ccdread',"rdnoise")
    iraf.datapars.setParam('exposure',"exptime")
    
    iraf.datapars.setParam('filter',"filter")
    iraf.datapars.setParam('obstime',"time-obs")
    iraf.datapars.setParam('sigma',"INDEF")
    iraf.photpars.setParam('zmag',0.)
    iraf.centerpars.setParam('cbox',9.)
    iraf.centerpars.setParam('maxshift',3.)
    iraf.fitskypars.setParam('salgorithm',"median")
    iraf.fitskypars.setParam('dannulus',10.)
    
    iraf.datapars.setParam('airmass','airmass')
    iraf.datapars.setParam('fwhmpsf',fwhm)
    iraf.photpars.setParam('apertures',3.*fwhm) # use a big aperture for this
    iraf.fitskypars.setParam('annulus',4.*fwhm)
    
    if not os.path.isfile(output):
        iraf.apphot.phot(image=img, coords=coords, output=output)
    with open(phot_tbl,'w+') as txdump_out :
        iraf.ptools.txdump(textfiles=output, fields="id,mag,merr,msky,stdev,rapert,xcen,ycen,ifilter,xairmass,peak,flux,image", expr='MAG != INDEF && MERR != INDEF', headers='no', Stdout=txdump_out)

    gMAG, gMERR, gSKY, gSERR, gRAPERT, gXPOS, gYPOS = np.loadtxt(phot_tbl, usecols=(1,2,3,4,5,6,7), dtype=float, unpack=True)
    gXAIRMASS = np.loadtxt(phot_tbl, usecols=(9,), dtype=float, unpack=True)
    gFILTER = np.loadtxt(phot_tbl, usecols=(8,), dtype=str, unpack=True)
    gID = np.loadtxt(phot_tbl, usecols=(0,), dtype=int, unpack=True)

    # keep the actual ID number to select from SDSS stars
    keep = gID - 1

    # read in the catalog values
    x, y, ra, dec, u, ue, g, ge, r, re, i, ie, z, ze = np.loadtxt(coords, usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13), unpack=True)

    # pick out the ones that match the good phot stars
    g, ge, r, re, i, ie = np.array(g[keep]), np.array(ge[keep]), np.array(r[keep]), np.array(re[keep]), np.array(i[keep]), np.array(ie[keep])

    # and reduce the other vectors
    # gXPOS, gYPOS, gMAG, gMERR, gSKY, gSERR = np.array(gXPOS[keep]), np.array(gYPOS[keep]), np.array(gMAG[keep]), np.array(gMERR[keep]), np.array(gSKY[keep]), np.array(gSERR[keep])

    gXAIRMASS = gXAIRMASS[0]
    gRAPERT = gRAPERT[0]

    podicut, sdsscut = 0.005, 0.02
    # pick the right airmass extinction coefficient
    kdict = {'odi_g':kg, 'odi_r':kr, 'odi_i':ki}
    mdict = {'odi_g':g, 'odi_r':r, 'odi_i':i}
    k = kdict[gFILTER[0]]
    mag_cat = mdict[gFILTER[0]]
    g0 = gMAG - k*gXAIRMASS
    zp = mag_cat - g0
    color = g - r
    errcut = [j for j in range(len(gMERR)) if (gMERR[j] < podicut and ge[j] < sdsscut)]

    print 'fitting wtih '+repr(len(g0[errcut]))+' stars...'

    # fit zero point
    # linear lsq with numpy.polyfit
    p, pcov = np.polyfit(color[errcut], zp[errcut], 1, cov=True)
    perr = np.sqrt(np.diag(pcov))
    eps_g, zp_g, std_eps_g, std_zp_g = p[0], p[1], perr[0], perr[1]

    print '--------------------------------------------------------------------------'
    print 'Here are the fit values:'
    print 'eps  '+'      std_eps  '+'  zp  '+'         std_zp'
    print '{0:10.7f} {1:10.7f} {2:10.7f} {3:10.7f}'.format(eps_g, std_eps_g, zp_g, std_zp_g)
    star_zp_g = zp - eps_g*color
    print 'std. dev. in ZP per star (not fit): {0:10.7f}'.format(np.std(star_zp_g[errcut]))
    # for i in range(len(zp)):
    #     print color[i], zp[i], gMERR[i], ge[i]

    print np.median(star_zp_g[errcut]), np.std(star_zp_g[errcut])
    return np.median(star_zp_g[errcut]), np.std(star_zp_g[errcut])

def main():
    pass

if __name__ == '__main__':
    main()

