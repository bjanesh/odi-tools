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

def trim_img(img):
    x1,x2 = 2508,15798
    y1,y2 = 2216,15506
    input = img[:-5]+'['+repr(x1)+':'+repr(x2)+','+repr(y1)+':'+repr(y2)+']'
    output =  img[:-5]+'.trim.fits'
    if not os.path.isfile(output):
	print 'Trimming image: ' ,img
        iraf.unlearn(iraf.imcopy)
        iraf.imcopy(input = input,output = output,verbose='no',mode='h')
    
     

def find_sources_full(img,fwhm,bg_std,threshold=4.0):
    output = img[:-10]+'_sources.coo'
    if not os.path.isfile(output):
        print 'Locating sources on ',img
        print 'Will output to ',output
        iraf.unlearn(iraf.apphot.daofind)
        iraf.datapars.setParam('fwhmpsf',fwhm,check=1)
        iraf.datapars.setParam('datamin',-900,check=1)
        iraf.datapars.setParam('datamax',60000,check=1)
        iraf.datapars.setParam('sigma',bg_std,check=1)
        iraf.findpars.setParam('threshold',threshold)
        iraf.apphot.daofind.setParam('output',output)
        iraf.apphot.daofind(image=img, verbose="no", verify='no')
        
def phot_sources_full(img,fwhm,airmass):
    iraf.ptools(_doprint=0)
    coords = img[:-10]+'_sources.coo'
    output = img[:-10]+'.phot.1'
    phot_tbl = img[0:-10]+'.srcphot'
    if not os.path.isfile(phot_tbl) :
        print 'phot-ing ', img, ' from daofind'
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
        iraf.photpars.setParam('apertures',2.*fwhm)
        iraf.fitskypars.setParam('annulus',4.*fwhm)

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
    phot_tbl = img[0:-10]+'.srcphot'
    outputradec = img[0:-10]+'.srcphotrd'
    hdulist= odi.fits.open(img)
    
    if inst == 'podi':
        pvlist = hdulist[0].header['PV*']
        for pv in pvlist:
            tpv = 'T'+pv
            hdulist[0].header.rename_keyword(pv, tpv, force=False)
    w = odi.WCS(hdulist[0].header)
    
    MAG, MERR, SKY, SERR, RAPERT, XPOS, YPOS = np.loadtxt(phot_tbl, usecols=(1,2,3,4,5,6,7), dtype=float, unpack=True)
    with open(outputradec, 'w+') as fxy:
        for i,c in enumerate(XPOS):
            coords2 = [[XPOS[i],YPOS[i]]]
            pixcrd2 = w.wcs_pix2world(coords2, 1)
            print >> fxy, pixcrd2[0][0], pixcrd2[0][1],XPOS[i],YPOS[i],MAG[i], MERR[i],SKY[i],SERR[i],RAPERT[i]
    hdulist.close()
    
def match_phot_srcs(img1,img2):

    img1_srcs =img1[0:-10]+'.srcphotrd'
    img1_srsc_match = img1[:-10]+'.match.srscrd'
    img2_srcs =img2[0:-10]+'.srcphotrd'
    img2_srsc_match = img2[:-10]+'.match.srscrd'
    
    ra_1, dec_1,x_1,y_1,mag_1,merr_1,sky_1,serr_1,rapert_1 = np.loadtxt(img1_srcs,usecols=(0,1,2,3,4,5,6,7,8),unpack=True)
    ra_2, dec_2,x_2,y_2,mag_2,merr_2,sky_2,serr_2,rapert_2 = np.loadtxt(img2_srcs,usecols=(0,1,2,3,4,5,6,7,8),unpack=True)
    

    img1_catalog = SkyCoord(ra = ra_1*u.degree, dec= dec_1*u.degree)
    
    img2_catalog = SkyCoord(ra = ra_2*u.degree, dec= dec_2*u.degree)
    
    id_img1, id_img2, d2d, d3d = img2_catalog.search_around_sky(img1_catalog,0.00001*u.deg)
    
    print len(id_img1),len(id_img2),len(x_1),len(x_2)
    
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
    
    junk = open('match.reg',"w")
    for i in range(len(ra_1)):
	print >> junk, x_1[i], y_1[i]
    junk.close()
    
    