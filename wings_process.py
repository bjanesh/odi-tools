#! /usr/local/bin/python
import os, sys, time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
# import matplotlib.pyplot as plt
from subprocess import call
# from astropy import wcs
# import sewpy
from astropy.io import fits
from pyraf import iraf
from escut_new import escut 
import odi_calibrate as odi
from matplotlib import gridspec

home_root = os.environ['HOME']

iraf.images(_doprint=0)
iraf.tv(_doprint=0)
iraf.ptools(_doprint=0)
iraf.noao(_doprint=0)
iraf.digiphot(_doprint=0)
iraf.photcal(_doprint=0)
iraf.apphot(_doprint=0)  
iraf.imutil(_doprint=0)

# set the title string
title_string = 'm13-se-short'

filter_id = np.loadtxt(title_string+'.aux', usecols=(1,), dtype=str, unpack=True)
exptime, airmass, fwhm_all_pix, fwhm_sdss_pix, skybgstd, skylevel = np.loadtxt(title_string+'.aux', usecols=(2,3,6,7,8,9), unpack=True)

r_exp = [i for i,f in enumerate(filter_id) if f.endswith('r')]
g_exp = [i for i,f in enumerate(filter_id) if f.endswith('g')]

# r_exp = np.where((exptime<150.0) & (exptime>50.0))
# g_exp = np.where(exptime>150.0)

airmass_g = np.mean(airmass[g_exp])
airmass_r = np.mean(airmass[r_exp])

fwhm_g = max(fwhm_sdss_pix[g_exp])
fwhm_r = max(fwhm_sdss_pix[r_exp])

bgm_g = np.mean(skylevel[g_exp])
bgm_r = np.mean(skylevel[r_exp])

bg_g = np.mean(skybgstd[g_exp])
bg_r = np.mean(skybgstd[r_exp])

odi.download_sdss(title_string+'.g.fits', title_string+'.r.fits', gmaglim = 21)
odi.calibrate(title_string+'.g.fits', title_string+'.r.fits', fwhm_g, fwhm_r, airmass_g, airmass_r)

if not os.path.isfile(title_string+'.r.sh.fits'):
    odi.imalign(title_string+'.g.fits', title_string+'.r.fits')

# copy the good region (no cell gaps) to a new file        
fits_g = title_string+'.g.sh.fits'
fits_r = title_string+'.r.sh.fits'
    
# make an imsets file
if not os.path.isfile(title_string+'.imsets') :
    imset_file = open(title_string+'.imsets', 'w+')
    print >> imset_file, title_string, ':', fits_r, fits_g
    imset_file.close()

kg = 0.200
kr = 0.12
ki = 0.058
rdnoise = 6.5

# get the photometric calibration coefficients from Steven's help file <--
# or from the image header/fits table/ whatever
photcalFile = open(title_string+'_help.txt')
photcal = photcalFile.read()
photcalLines = photcal.splitlines()

mu_gr = float(photcalLines[28].split()[5])
zp_gr = float(photcalLines[30].split()[4])
eps_gr = float(photcalLines[34].split()[5])
zp_r = float(photcalLines[36].split()[4])
amg = float(photcalLines[25].split()[5])
amr = float(photcalLines[26].split()[5])
photcalFile.close()

# print mu_gi, zp_gi, eps_gi, zp_r, amg, ami

fits_h_r = fits.open(fits_r)
fits_h_g = fits.open(fits_g)

# fwhm_i = fits_h_i[0].header['F_AVGSEE']/0.11
# fwhm_g = fits_h_g[0].header['F_AVGSEE']/0.11

# get steven's/QR's estimate of the image FWHMPSF
# try:
#     fwhm_i = fits_h_i[0].header['FWHMPSF']
#     fwhm_g = fits_h_g[0].header['FWHMPSF']
# except:
#     fwhm_i = fits_h_i[0].header['SEEING']/0.11
#     fwhm_g = fits_h_g[0].header['SEEING']/0.11

# for now, assume that the FWHMPSF is 10 pixels until we remeasure
# fwhm_g, fwhm_r = 10.0, 10.0

print 'Target Coordinates :: ',fits_h_r[0].header['TARGRA'],fits_h_r[0].header['TARGDEC']
print 'Image header FWHM :: g = {0:5.3f} : i = {1:5.3f}'.format(fwhm_g,fwhm_r)


print 'Image mean BG sigma value :: g = {0:5.3f} : i = {1:5.3f}'.format(bg_g,bg_r)
print 'Image mean BG median value :: g = {0:5.3f} : i = {1:5.3f}'.format(bgm_g,bgm_r)

# daofind steps
# find all the sources in the image (threshold value will be data dependent, 4.0 is good for UCHVCs)
# g image
if not os.path.isfile(fits_g+'.coo.1') :
    iraf.datapars.setParam('fwhmpsf',fwhm_g,check=1)
    iraf.datapars.setParam('sigma',bg_g,check=1)
    
    iraf.findpars.setParam('threshold',4.0)
    iraf.apphot.daofind(image=fits_g, verbose="no", verify='no')
#     
#     # i image
if not os.path.isfile(fits_r+'.coo.1') :
    iraf.datapars.setParam('fwhmpsf',fwhm_r,check=1)
    iraf.datapars.setParam('sigma',bg_r,check=1)
    
    iraf.findpars.setParam('threshold',4.0)
    iraf.apphot.daofind(image=fits_r, verbose="no", verify='no')

# now phot the stars found in daofind
if not os.path.isfile(fits_g+'.mag.1') :
    print 'phot-ing g band daofind stars. This is going to take a while...'
    iraf.unlearn(iraf.apphot.phot,iraf.datapars,iraf.photpars,iraf.centerpars,iraf.fitskypars)
    iraf.apphot.phot.setParam('interactive',"no")
    iraf.apphot.phot.setParam('verify',"no")
    iraf.datapars.setParam('datamax',50000.)
    iraf.datapars.setParam('gain',"gain")
    iraf.datapars.setParam('ccdread',rdnoise)
    iraf.datapars.setParam('exposure',"exptime")
    iraf.datapars.setParam('airmass',airmass_g)
    iraf.datapars.setParam('filter',"filter")
    iraf.datapars.setParam('obstime',"time-obs")
    iraf.datapars.setParam('sigma',"INDEF")
    iraf.photpars.setParam('zmag',0.)
    iraf.centerpars.setParam('cbox',9.)
    iraf.centerpars.setParam('maxshift',3.)
    iraf.fitskypars.setParam('salgorithm',"median")
    iraf.fitskypars.setParam('dannulus',10.)
    
    iraf.datapars.setParam('fwhmpsf',fwhm_g)
    iraf.photpars.setParam('apertures',2.*fwhm_g)
    iraf.fitskypars.setParam('annulus',4.*fwhm_g)

    iraf.apphot.phot(image=fits_g, coords=fits_g+'.coo.1')

if not os.path.isfile(fits_r+'.mag.1') :
    print 'phot-ing r band daofind stars. This is going to take a while...'
    iraf.unlearn(iraf.phot,iraf.datapars,iraf.photpars,iraf.centerpars,iraf.fitskypars)
    iraf.apphot.phot.setParam('interactive',"no")
    iraf.apphot.phot.setParam('verify',"no")
    iraf.datapars.setParam('datamax',50000.)
    iraf.datapars.setParam('gain',"gain")
    iraf.datapars.setParam('ccdread',rdnoise)
    iraf.datapars.setParam('exposure',"exptime")
    iraf.datapars.setParam('airmass',airmass_r)
    iraf.datapars.setParam('filter',"filter")
    iraf.datapars.setParam('obstime',"time-obs")
    iraf.datapars.setParam('sigma',"INDEF")
    iraf.photpars.setParam('zmag',0.)
    iraf.centerpars.setParam('cbox',9.)
    iraf.centerpars.setParam('maxshift',3.)
    iraf.fitskypars.setParam('salgorithm',"median")
    iraf.fitskypars.setParam('dannulus',10.)
    iraf.datapars.setParam('fwhmpsf',fwhm_r)
    iraf.photpars.setParam('apertures',2.*fwhm_r)
    iraf.fitskypars.setParam('annulus',4.*fwhm_r)

    iraf.apphot.phot(image=fits_r, coords=fits_r+'.coo.1')

# get rid of regions you don't want using pselect    
# if not os.path.isfile('mask.reg'):
#     print 'To continue you should mask out bright stars, galaxies, etc.'
#     print 'in DS9 and export to an IRAF PROS file named mask.reg'
#     raw_input("Press Enter when finished:")


# m3,m4,m5,m6 = np.loadtxt('mask.reg',usecols=(2,3,4,5),unpack=True)
# # g image pselects
if not os.path.isfile(fits_g+'.mag.1a') :
    if os.path.isfile('temp2') :    
        os.remove('temp2')
    iraf.ptools.pselect(infi=fits_g+'.mag.1', outfi='temp1', expr="MAG != INDEF")
#     for i in range(len(m3)) :
#         mx1 = m3[i] - (m5[i]/2.)
#         mx2 = m3[i] + (m5[i]/2.)
#         my1 = m4[i] - (m6[i]/2.)
#         my2 = m4[i] + (m6[i]/2.)
#         iraf.ptools.pselect(infi='temp1', outfi='temp2', expr='(XCE < '+repr(int(mx1))+' || XCE > '+repr(int(mx2))+') || (YCE < '+repr(int(my1))+' || YCE > '+repr(int(my2))+')')
#         os.rename('temp2', 'temp1')
    os.rename('temp1', fits_g+'.mag.1a')
    if os.path.isfile('temp1') :    
        os.remove('temp1')
# 
# # i image pselects
if not os.path.isfile(fits_r+'.mag.1a') :
    if os.path.isfile('temp2') :    
        os.remove('temp2')
    iraf.ptools.pselect(infi=fits_r+'.mag.1', outfi='temp1', expr="MAG != INDEF")
#     for i in range(len(m3)) :
#         mx1 = m3[i] - (m5[i]/2.)
#         mx2 = m3[i]+ (m5[i]/2.)
#         my1 = m4[i] - (m6[i]/2.)
#         my2 = m4[i] + (m6[i]/2.)
#         iraf.ptools.pselect(infi='temp1', outfi='temp2', expr='(XCE < '+repr(int(mx1))+' || XCE > '+repr(int(mx2))+') || (YCE < '+repr(int(my1))+' || YCE > '+repr(int(my2))+')')
#         os.rename('temp2', 'temp1')
    os.rename('temp1', fits_r+'.mag.1a')
    if os.path.isfile('temp1') :    
        os.remove('temp1')

# mkobsfile stuff
# mkobsfile MATCHES sources between images to get rid of random sources, things that are masked in one or the other, etc.
iraf.digiphot.mkobsfile.setParam('photfiles',title_string+'.*.fits.mag.1a')
iraf.digiphot.mkobsfile.setParam('idfilters','odi_r,odi_g')
iraf.digiphot.mkobsfile.setParam('imsets',title_string+'.imsets')
iraf.digiphot.mkobsfile.setParam('obscolumns','2 3 4 5')
iraf.digiphot.mkobsfile.setParam('shifts',None)
iraf.digiphot.mkobsfile.setParam('apercors',None)
iraf.digiphot.mkobsfile.setParam('allfilters','yes')

# if not os.path.isfile('ifirst_tol6.out') :
#     iraf.mkobsfile.setParam('observations','ifirst_tol6.out')
#     iraf.mkobsfile.setParam('tolerance',6.)
#     iraf.mkobsfile()

if not os.path.isfile('mkobsfile.out') :
    iraf.digiphot.mkobsfile.setParam('observations','mkobsfile.out')
    iraf.digiphot.mkobsfile.setParam('tolerance',7.) # number of pixels away matched source can be, DATA DEPENDENT!
    iraf.digiphot.mkobsfile()

# if not os.path.isfile('ifirst_tol8.out') :
#     iraf.mkobsfile.setParam('observations','ifirst_tol8.out')
#     iraf.mkobsfile.setParam('tolerance',8.)
#     iraf.mkobsfile()
# call("awk '{ if ($2 ~ "odi_r") print $5, $6 }' ifirst_tol7.dat > tol7_r.pos")

# print matched sources to a file suitable for marking
if os.path.isfile('mkobsfile.out') :
    mx,my = np.loadtxt('mkobsfile.out',usecols=(4,5),unpack=True)
    mfilter = np.loadtxt('mkobsfile.out',usecols=(1,),dtype=str,unpack=True)
    match_pos_file_g = open("match_g.pos", 'w+')
    match_pos_file_r = open("match_r.pos", 'w+')
    for i in range(len(mx)) :
        if mfilter[i]== 'odi_g' :
            print >> match_pos_file_g, mx[i], my[i]
        if mfilter[i] == 'odi_r' :
            print >> match_pos_file_r, mx[i], my[i]
    match_pos_file_g.close()
    match_pos_file_r.close()
    
def getfwhm(image, coords, outputfile, radius=4.0, buff=7.0, width=5.0, rplot=15.0, center='yes'):
    '''
    Get a fwhm estimate for the image using the SDSS catalog stars and IRAF imexam (SLOW, but works)
    Adapted from Kathy's getfwhm script (this implementation is simpler in practice)
    '''
    iraf.tv.rimexam.setParam('radius',radius)
    iraf.tv.rimexam.setParam('buffer',buff)
    iraf.tv.rimexam.setParam('width',width)
    iraf.tv.rimexam.setParam('rplot',rplot)
    iraf.tv.rimexam.setParam('center',center)
    # fit a gaussian, rather than a moffat profile (it's more robust for faint sources)
    iraf.tv.rimexam.setParam('fittype','gaussian')
    iraf.tv.rimexam.setParam('iterati',1)

    if not os.path.isfile(outputfile):
        iraf.tv.imexamine( image, logfile = outputfile, keeplog = 'yes', defkey = "a", imagecur = coords, wcs = "world", use_display='no')

# you might want to remeasure the FWHMs to get a better global estimate now that we (should) only have good sources in the image
# use getfwhm, which is just a loop on imexam. (try to improve this with ralf's qr code)
if not os.path.isfile('getfwhm_g.log') :
    getfwhm(fits_g, 'match_g.pos', 'getfwhm_g.log')
# 
if not os.path.isfile('getfwhm_r.log') :
    getfwhm(fits_r, 'match_r.pos', 'getfwhm_r.log')
    
# determine the aperture correction needed--this is actually an extremely important step. uses 4.5x the measured FWHM as the aperture DATA DEPENDENT 
if not os.path.isfile('apcor.tbl.txt'):
    ap_gx,ap_gy = np.loadtxt('getfwhm_g.log',usecols=(0,1),unpack=True)
    ap_mag_g = np.loadtxt('getfwhm_g.log',usecols=(5,),dtype=str,unpack=True)
    ap_peak_g = np.loadtxt('getfwhm_g.log',usecols=(8,),dtype=str,unpack=True)
    ap_fwhm_g = np.loadtxt('getfwhm_g.log',usecols=(12,),dtype=str,unpack=True)
    ap_rx,ap_ry = np.loadtxt('getfwhm_r.log',usecols=(0,1),unpack=True)
    ap_mag_r = np.loadtxt('getfwhm_r.log',usecols=(5,),dtype=str,unpack=True)
    ap_peak_r = np.loadtxt('getfwhm_r.log',usecols=(8,),dtype=str,unpack=True)
    ap_fwhm_r = np.loadtxt('getfwhm_r.log',usecols=(12,),dtype=str,unpack=True)
    
    ap_cand1_g = [(ap_gx[i],ap_gy[i],float(ap_fwhm_g[i]),float(ap_peak_g[i]),float(ap_mag_g[i])) for i in range(len(ap_gx)) if (ap_peak_g[i] != 'INDEF' and ap_fwhm_g[i] != 'INDEF' and ap_mag_g[i] != 'INDEF')]
    ap_cand1_r = [(ap_rx[i],ap_ry[i],float(ap_fwhm_r[i]),float(ap_peak_r[i]),float(ap_mag_r[i])) for i in range(len(ap_rx)) if (ap_peak_r[i] != 'INDEF' and ap_fwhm_r[i] != 'INDEF' and ap_mag_r[i] != 'INDEF')]
    
    if fwhm_r < 11.0 :
        ap_cand_g = [ap_cand1_g[i] for i in range(len(ap_cand1_g)) if (20000. < ap_cand1_g[i][3] < 50000.)]
        ap_cand_r = [ap_cand1_r[i] for i in range(len(ap_cand1_r)) if (20000. < ap_cand1_r[i][3] < 50000.)]
        # print ap_cand_g
        ap_avg_g1 = np.mean([ap_cand_g[i][2] for i in range(len(ap_cand_g))])
        ap_avg_r1 = np.mean([ap_cand_r[i][2] for i in range(len(ap_cand_r))])
    
        ap_std_g1 = np.std([ap_cand_g[i][2] for i in range(len(ap_cand_g))])
        ap_std_r1 = np.std([ap_cand_r[i][2] for i in range(len(ap_cand_r))])
        # print ap_avg_g1,ap_avg_r1,ap_std_g1,ap_std_r1
        ap_stars_g = [ap_cand_g[i] for i in range(len(ap_cand_g)) if ((ap_avg_g1-ap_std_g1) < ap_cand_g[i][2] < (ap_avg_g1+ap_std_g1))]
        ap_stars_r = [ap_cand_r[i] for i in range(len(ap_cand_r)) if ((ap_avg_r1-ap_std_r1) < ap_cand_r[i][2] < (ap_avg_r1+ap_std_r1))]
        # print len(ap_stars_g), len(ap_stars_r)
        ap_avg_g = np.mean([ap_stars_g[i][2] for i in range(len(ap_stars_g))])
        ap_avg_r = np.mean([ap_stars_r[i][2] for i in range(len(ap_stars_r))])
    
        ap_std_g = np.std([ap_stars_g[i][2] for i in range(len(ap_stars_g))])
        ap_std_r = np.std([ap_stars_r[i][2] for i in range(len(ap_stars_r))])
        print 'Measured image FWHM :: g = {0:5.3f} : i = {1:5.3f}'.format(ap_avg_g,ap_avg_r)
    
        if not os.path.isfile('apcor_stars_g.txt'):
            ap_file_g = open('apcor_stars_g.txt','w+')
            for i in range(len(ap_stars_g)) :
                print >> ap_file_g, ap_stars_g[i][0], ap_stars_g[i][1], ap_stars_g[i][2], ap_stars_g[i][3]
            ap_file_g.close()
        
        if not os.path.isfile('apcor_stars_r.txt'):
            ap_file_r = open('apcor_stars_r.txt','w+')
            for i in range(len(ap_stars_r)) :
                print >> ap_file_r, ap_stars_r[i][0], ap_stars_r[i][1], ap_stars_r[i][2], ap_stars_r[i][3]
            ap_file_r.close()
    
        iraf.unlearn(iraf.apphot.phot,iraf.datapars,iraf.photpars,iraf.centerpars,iraf.fitskypars)
    
        iraf.apphot.phot.setParam('interactive','no')
        iraf.apphot.phot.setParam('verify','no')
    
        iraf.datapars.setParam('gain',"gain")
        iraf.datapars.setParam('ccdread',"rdnoise")
        iraf.datapars.setParam('exposure',"exptime")
        iraf.datapars.setParam('airmass',"airmass")
        iraf.datapars.setParam('filter',"filter")
        iraf.datapars.setParam('obstime',"time-obs")
        iraf.datapars.setParam('sigma','INDEF')
        iraf.photpars.setParam('zmag',0.)
        iraf.centerpars.setParam('cbox',9.)
        iraf.centerpars.setParam('maxshift',3.)
        iraf.fitskypars.setParam('salgorithm',"median")
        iraf.fitskypars.setParam('dannulus',10.)
    
        iraf.datapars.setParam('fwhmpsf',ap_avg_g)
        iraf.photpars.setParam('apertures',str('"'+repr(ap_avg_g)+','+repr(5.0*ap_avg_g)+'"'))
        iraf.fitskypars.setParam('annulus',6.5*ap_avg_g)
    
        iraf.apphot.phot(image=fits_g, coords='apcor_stars_g.txt', output="g.apcor.mag.1")
    
        iraf.datapars.setParam('fwhmpsf',ap_avg_r)
        iraf.photpars.setParam('apertures',str('"'+repr(ap_avg_r)+','+repr(5.0*ap_avg_r)+'"'))
        iraf.fitskypars.setParam('annulus',6.5*ap_avg_r)
        iraf.apphot.phot(image=fits_r, coords='apcor_stars_r.txt', output="i.apcor.mag.1")
    
        apt_file_g = open("apcor_table_g.txt", 'w+')
        apt_file_r = open("apcor_table_r.txt", 'w+')
    
        iraf.ptools.txdump(textfiles='g.apcor.mag.1', fields="ID,XCEN,YCEN,MAG", expr='yes', Stdout=apt_file_g)
        iraf.ptools.txdump(textfiles='i.apcor.mag.1', fields="ID,XCEN,YCEN,MAG", expr='yes', Stdout=apt_file_r)
        apt_file_g.close()
        apt_file_r.close()
    
        onex_g,four5x_g = np.loadtxt('apcor_table_g.txt',usecols=(3,4),unpack=True)
        onex_r,four5x_r = np.loadtxt('apcor_table_r.txt',usecols=(3,4),unpack=True)
    
        apcor_rnd_g = four5x_g - onex_g
        apcor_g = np.mean(apcor_rnd_g)
        apcor_std_g = np.std(apcor_rnd_g)
        apcor_sem_g = apcor_std_g/np.sqrt(len(apcor_rnd_g))
    
        apcor_rnd_r = four5x_r - onex_r
        apcor_r = np.mean(apcor_rnd_r)
        apcor_std_r = np.std(apcor_rnd_r)
        apcor_sem_r = apcor_std_r/np.sqrt(len(apcor_rnd_r))
    
        print 'Aperture correction :: g = {0:7.4f} : i = {1:7.4f}'.format(apcor_g,apcor_r)
        print 'Aperture corr. StD. :: g = {0:6.4f} : i = {1:6.4f}'.format(apcor_std_g,apcor_std_r)
        print 'Aperture corr. SEM  :: g = {0:6.4f} : i = {1:6.4f}'.format(apcor_sem_g,apcor_sem_r)
        print 'Aperture corr. N    :: g = {0:2d} : i = {1:2d}'.format(len(onex_g),len(onex_r))
    
        apcor_tbl = open('apcor.tbl.txt','w+')
        print >> apcor_tbl, apcor_g, apcor_std_g, apcor_sem_g
        print >> apcor_tbl, apcor_r, apcor_std_r, apcor_sem_r
        apcor_tbl.close()
    else :
        apcor_g = 0.0
        apcor_r = 0.0
        print 'Seeing is pretty bad, no aperture correction applied.'
else :
    apcor, apcor_std, apcor_sem = np.loadtxt('apcor.tbl.txt', usecols=(0,1,2), unpack=True)
    apcor_g = apcor[0]
    apcor_r = apcor[1]
    apcor_std_g = apcor_std[0]
    apcor_std_r = apcor_std[1]
    apcor_sem_g = apcor_sem[0]
    apcor_sem_r = apcor_sem[1]
    print 'Aperture correction :: g = {0:7.4f} : i = {1:7.4f}'.format(apcor_g,apcor_r)
    print 'Aperture corr. StD. :: g = {0:6.4f} : i = {1:6.4f}'.format(apcor_std_g,apcor_std_r)
    print 'Aperture corr. SEM  :: g = {0:6.4f} : i = {1:6.4f}'.format(apcor_sem_g,apcor_sem_r)

    ap_gx,ap_gy = np.loadtxt('getfwhm_g.log',usecols=(0,1),unpack=True)
    ap_mag_g = np.loadtxt('getfwhm_g.log',usecols=(5,),dtype=str,unpack=True)
    ap_peak_g = np.loadtxt('getfwhm_g.log',usecols=(8,),dtype=str,unpack=True)
    ap_fwhm_g = np.loadtxt('getfwhm_g.log',usecols=(12,),dtype=str,unpack=True)
    good = np.where(ap_fwhm_g!='INDEF')
    fwhm_good = ap_fwhm_g[good].astype(float)
    ap_avg_g = np.median(fwhm_good)

    ap_rx,ap_ry = np.loadtxt('getfwhm_r.log',usecols=(0,1),unpack=True)
    ap_mag_r = np.loadtxt('getfwhm_r.log',usecols=(5,),dtype=str,unpack=True)
    ap_peak_r = np.loadtxt('getfwhm_r.log',usecols=(8,),dtype=str,unpack=True)
    ap_fwhm_r = np.loadtxt('getfwhm_r.log',usecols=(12,),dtype=str,unpack=True)
    good = np.where(ap_fwhm_r!='INDEF')
    fwhm_good = ap_fwhm_r[good].astype(float)
    ap_avg_r = np.median(fwhm_good)
    
# try :
#     if not os.path.isfile('point_source_rnfo') :
#         call('escut')
# except :
#     print 'Something went wrong running escut. Try running it outside of this script.'

# escut_g = escut(fits_g, 'match_g.pos', ap_fwhm_g)
# print ap_peak_r.size, ap_fwhm_r.size
escut_r = escut(fits_r, 'match_r.pos', ap_fwhm_r, ap_peak_r)
    
# finally rephot just the good stuff to get a good number
# if os.path.isfile('escut_g.pos') :
iraf.unlearn(iraf.phot,iraf.datapars,iraf.photpars,iraf.centerpars,iraf.fitskypars)
#
iraf.apphot.phot.setParam('interactive',"no")
iraf.apphot.phot.setParam('verify',"no")
iraf.datapars.setParam('datamin',"INDEF")
iraf.datapars.setParam('datamax',50000.)
iraf.datapars.setParam('gain',"gain")
iraf.datapars.setParam('ccdread',"rdnoise")
iraf.datapars.setParam('exposure',"exptime")
iraf.datapars.setParam('airmass',"airmass")
iraf.datapars.setParam('filter',"filter")
iraf.datapars.setParam('obstime',"time-obs")
iraf.datapars.setParam('sigma',"INDEF")
iraf.photpars.setParam('zmag',0.)
iraf.centerpars.setParam('calgorithm',"centroid")
iraf.centerpars.setParam('cbox',9.)
iraf.centerpars.setParam('maxshift',3.)
iraf.fitskypars.setParam('salgorithm',"median")
iraf.fitskypars.setParam('dannulus',10.)
#
# Use an aperture that is 1 x <fwhm>, because an aperture correction
# will be applied in the calc_calib_mags step
# Using a sky annulus thatbegins at 6 x <fwhm> should be fine
# g-band
if not os.path.isfile(title_string+'_sources_g.mag.1') :
    print 'Phot-ing g band point sources, this could take a while.'
    iraf.datapars.setParam('fwhmpsf',ap_avg_g)
    iraf.photpars.setParam('apertures',ap_avg_g)
    iraf.fitskypars.setParam('annulus',6*ap_avg_g)
    iraf.apphot.phot(image=fits_g, coords='escut_g.pos', output=title_string+'_sources_g.mag.1')

# i-band    
# if os.path.isfile('escut_r.pos') :
if not os.path.isfile(title_string+'_sources_r.mag.1') :
    print 'Phot-ing r band point sources, this could take a while.'
    iraf.datapars.setParam('fwhmpsf',ap_avg_r)
    iraf.photpars.setParam('apertures',ap_avg_r)
    iraf.fitskypars.setParam('annulus',6*ap_avg_r)
    iraf.apphot.phot(image=fits_r, coords='escut_r.pos', output=title_string+'_sources_r.mag.1')


def getexttbl(ra,dec,fname='extinction.tbl.txt'):
    '''Takes RA and Dec in colon separated sexagesimal format and queries 
        http://irsa.ipac.caltech.edu/ for the extinction table at that location
        File is saved as extinction.tbl.txt or as specified by user in optional 3rd arg'''
    import urllib2
    from bs4 import BeautifulSoup
    # parse the ra and dec
    ravals = ra.split(':')
    decvals = dec.split(':')

    ra_str = ravals[0]+'h'+ravals[1]+'m'+ravals[2]+'s'
    dec_str = decvals[0]+'d'+decvals[1]+'m'+decvals[2]+'s'

    # set the url to go get the extinction table
    exturl = "http://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr="+ra_str+'+'+dec_str+'+equ+j2000' 
    print '-> from', exturl

    xmlget = urllib2.urlopen(exturl)
    soup = BeautifulSoup(xmlget,"lxml-xml")
    # print soup.prettify()

    tblurl = soup.result.data.table.string
    exttbl = urllib2.urlopen(tblurl)

    f = open(fname,'w+')
    print >> f, exttbl.read()
    f.close()

if not os.path.isfile('extinction.tbl.txt'):
    print 'Fetching extinction table for',fits_h_r[0].header['TARGRA'],fits_h_r[0].header['TARGDEC']
    getexttbl(fits_h_r[0].header['TARGRA'],fits_h_r[0].header['TARGDEC'])

LamEff,A_over_E_B_V_SandF,A_SandF,A_over_E_B_V_SFD,A_SFD= np.genfromtxt('extinction.tbl.txt', usecols=(2,3,4,5,6),unpack=True,skip_header=27,skip_footer=12)
A_rd = np.genfromtxt('extinction.tbl.txt', usecols=(1,),dtype=str,unpack=True,skip_header=27,skip_footer=12)
E_B_V = np.genfromtxt('extinction.tbl.txt', usecols=(2,),skip_header=1,skip_footer=42)

for j in range(len(A_rd)):
    if A_rd[j] == 'g':
        cal_A_g = A_over_E_B_V_SandF[j]*0.86*E_B_V # E(B-V) is the Schlegel+ value, S&F say with their calibration
for j in range(len(A_rd)):                                  # use 0.86*E(B-V) instead. cf. S&F2011 pg 1, 2011ApJ...737..103S
    if A_rd[j] == 'i':
        cal_A_r = A_over_E_B_V_SandF[j]*0.86*E_B_V
        
print 'Reddening correction :: g = {0:7.4f} : i = {1:7.4f}'.format(cal_A_g,cal_A_r)

txdump_out = open('phot_sources.txdump','w+')
iraf.ptools.txdump(textfiles=title_string+'_sources_*.mag.1', fields="id,mag,merr,msky,stdev,rapert,xcen,ycen,ifilter,xairmass,image", expr='yes', headers='no', Stdout=txdump_out)
txdump_out.close()

call('sort -g phot_sources.txdump > temp', shell=True)
call('mv temp phot_sources.txdump', shell=True)
call('awk -f ~/projects/wings/make_calibdat phot_sources.txdump > calibration.dat', shell=True)

nid,gx,gy,g_i,g_ierr,rx,ry,r_i,r_ierr = np.loadtxt('calibration.dat',usecols=(0,1,2,4,5,6,7,9,10),unpack=True)

g0 = g_i - (kg*amg) + apcor_g 
r0 = r_i - (kr*amr) + apcor_r
gmr = mu_gr*(g0-r0) + zp_gr

r_mag = r0 + eps_gr*gmr + zp_r
g_mag = gmr + r_mag - cal_A_g 
r_mag = r_mag - cal_A_r
gmr = g_mag - r_mag

print 'Median (g-r) :: g - r = {0:7.4f}'.format(np.median(gmr))
print 'Final number of phot-ed stars :: g = {0:5d} : r = {1:5d}'.format(len(g_mag),len(r_mag))

g_mag_lims = [g_mag[i] for i in range(len(g_mag)) if (g_ierr[i] >= 0.2)]
r_mag_lims = [r_mag[i] for i in range(len(r_mag)) if (r_ierr[i] >= 0.2)]

# print '5-sigma limit :: g = {0:7.4f} : i = {1:7.4f}'.format(min(g_mag_lims), min(i_mag_lims))

f3 = open('calibrated_mags.dat', 'w+')
for i in range(len(rx)) :
    print >> f3, '{0:8.2f} {1:8.2f} {2:12.3f} {3:12.3f} {4:8.2f} {5:8.2f} {6:12.3f} {7:12.3f} {8:12.3f} '.format(gx[i],gy[i],g_mag[i],g_ierr[i],rx[i],ry[i],r_mag[i],r_ierr[i],gmr[i])
f3.close()

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
plt.savefig(title_string+"_CMD.pdf")

# delete the pipeline WCS keywords from the header, Steven's are better
# iraf.imutil.hedit(images=fits_g, fields='PV*', delete='yes', verify='no')
# iraf.imutil.hedit(images=fits_r, fields='PV*', delete='yes', verify='no')

