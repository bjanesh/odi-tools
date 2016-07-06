import sys, os, glob, string
import pandas as pd
from collections import OrderedDict
import shutil
import numpy as np
from pyraf import iraf
from photutils.detection import detect_sources
from photutils import source_properties, properties_table
import odi_config as odi
import pandas as pd
import astropy
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting

def source_find(img,ota,inst):
    """
    This function will find sources on an OTA
    using the detect_sources module from photutils.
    This will return of csv file of the sources found
    with the x,y,Ra,Dec,source_sum,max_value, and
    elongation of the source. The elongation parameter is
    semimajor_axis / semiminor_axis. This output is needed
    for the source_xy function.
    """
    image = odi.reprojpath+'reproj_'+ota+'.'+str(img[16:])
    QR_raw = odi.fits.open(image)
    hdu_ota = QR_raw[0]

    if inst == 'podi':
        pvlist = hdu_ota.header['PV*']
        for pv in pvlist:
            tpv = 'T'+pv
            hdu_ota.header.rename_keyword(pv, tpv, force=False)
    w = odi.WCS(hdu_ota.header)
    #w.wcs.ctype = ["RA---TPV", "DEC--TPV"]
    bg_mean,bg_median,bg_std = odi.mask_ota(img,ota,reproj=True)
    threshold = bg_median + (bg_std * 5.)
    print bg_mean,bg_median,bg_std
    segm_img = detect_sources(hdu_ota.data, threshold, npixels=20)
    source_props = source_properties(hdu_ota.data,segm_img,wcs=w)

    columns = ['id', 'xcentroid', 'ycentroid', 'ra_icrs_centroid',
	       'dec_icrs_centroid','source_sum','max_value','elongation']
    source_tbl = properties_table(source_props,columns=columns)
    source_tbl_df = source_tbl.to_pandas()

    outputfile = odi.sourcepath+'source_'+ota+'.'+str(img[16:-5])+'.csv'

    source_tbl_df.to_csv(outputfile,index=False)
    QR_raw.close()

def source_xy(img,ota,gapmask,filter,inst):
    """
    This function will return the x,y positions of
    sources found by source_find that are not too
    close to gaps or the edges of the ota.
    """
    image = odi.reprojpath+'reproj_'+ota+'.'+str(img[16:])
    #image = odi.bgsubpath+'bgsub_'+ota+'.'+str(img[16:])
    input_xy = odi.sourcepath+'source_'+ota+'.'+str(img[16:-5])+'.csv'
    outputxy = odi.sourcepath+'source_'+ota+'.'+str(img[16:-5])+'.xy'
    id,xcentroid,ycentroid,ra_icrs_centroid,dec_icrs_centroid,source_sum,max_value,elongation = np.loadtxt(input_xy,usecols=(0,1,2,3,4,5,6,7), unpack=True, delimiter=',', skiprows=1)
    QR_raw = odi.fits.open(image)
    hdu_ota = QR_raw[0]

    if inst == 'podi':
        pvlist = hdu_ota.header['PV*']
        for pv in pvlist:
            tpv = 'T'+pv
            hdu_ota.header.rename_keyword(pv, tpv, force=False)

    w = odi.WCS(hdu_ota.header)
    xdim = hdu_ota.header['NAXIS1']
    ydim = hdu_ota.header['NAXIS2']

    with open(outputxy, 'w+') as fxy:
        for i,c in enumerate(xcentroid):
            coords2 = [[xcentroid[i],ycentroid[i]]]
            pixcrd2 = coords2
            if  100.0 <= pixcrd2[0][0] < xdim-100.0 and 100.0 <= pixcrd2[0][1] < ydim-100.0  and elongation[i] <=1.75:
                # make an image cutout of the gap mask
                x, y = int(round(pixcrd2[0][0])), int(round(pixcrd2[0][1]))
                cutout = gapmask[y-30:y+30,x-30:x+30]
                if not (cutout.astype(bool)).any():
                    print >> fxy, pixcrd2[0][0], pixcrd2[0][1], id[i],ra_icrs_centroid[i],dec_icrs_centroid[i],source_sum[i],max_value[i],elongation[i]
    QR_raw.close()
    fxy.close()

def getfwhm_source(img, ota, radius=4.0, buff=7.0, width=5.0):
    """
    Measure the FWHM using IRAF at the x,y positions that are
    returned by source_xy.
    """
    image = odi.reprojpath+'reproj_'+ota+'.'+str(img[16:])
    coords = odi.sourcepath+'source_'+ota+'.'+str(img[16:-5])+'.xy'
    print image, coords
    outputfile = odi.sourcepath+img[0:-5]+'.'+ota+'.fwhm.log'

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
    sfwhm = np.median(gfwhm[np.where(gfwhm < 900.0)])

    print 'median gwfhm in ota',ota+': ',sfwhm,'pixels'# (determined via QR)'
    return sfwhm

def phot_sources(img, ota, fwhm):
    """
    Run IRAF phot on the the sources found.
    """
    iraf.ptools(_doprint=0)
    # values determined by ralf/daniel @ wiyn
    kg = 0.20
    kr = 0.12
    ki = 0.058

    image = odi.reprojpath+'reproj_'+ota+'.'+str(img[16:])
    coords = odi.sourcepath+'source_'+ota+'.'+str(img[16:-5])+'.xy'
    output = odi.sourcepath+img[0:-5]+'.'+ota+'.phot.1'
    phot_tbl = odi.sourcepath+img[0:-5]+'.'+ota+'.sourcephot'

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
    outputfile_clean = open(phot_tbl.replace('.sourcephot','_clean.sourcephot'),"w")
    for line in open(phot_tbl,"r"):
        if not 'INDEF' in line:
            outputfile_clean.write(line)
        if 'INDEF' in line:
            outputfile_clean.write(line.replace('INDEF','999'))
    outputfile_clean.close()
    os.rename(phot_tbl.replace('.sourcephot','_clean.sourcephot'),phot_tbl)
    return phot_tbl

def phot_combine(img, ota):
    """
    Combine all of the information gather on the found sources.
    These will be all of the values returned by source_find,
    source_xy, getfwhm_source, phot_sources.
    """
    coords = odi.sourcepath+'source_'+ota+'.'+str(img[16:-5])+'.xy'

    x, y, id,ra_icrs_centroid,dec_icrs_centroid,source_sum,max_value,elongation = np.loadtxt(coords,usecols=(0,1,2,3,4,5,6,7),unpack=True)

    phot_tbl = odi.sourcepath+img[0:-5]+'.'+ota+'.sourcephot'

    MAG, MERR, SKY, SERR, RAPERT, XPOS, YPOS = np.loadtxt(phot_tbl, usecols=(1,2,3,4,5,6,7), dtype=float, unpack=True)

    fwhmfile = odi.sourcepath+img[0:-5]+'.'+ota+'.fwhm.log'

    peak,fwhm = np.loadtxt(fwhmfile, usecols=(9,10), unpack=True)

    output = odi.sourcepath+img[0:-5]+'.'+ota+'.totphot'

    with open(output, 'w+') as xy:
	for i in range(len(x)):
	    print >> xy,x[i], y[i], id[i],ra_icrs_centroid[i],dec_icrs_centroid[i],source_sum[i],max_value[i],elongation[i], MAG[i], MERR[i], SKY[i], SERR[i], RAPERT[i], XPOS[i], YPOS[i], fwhm[i],peak[i]
    xy.close()

def source_scale(img,ref,filter):
    """
    This function calculates the scaling based on a reference image.
    The tables returned by phot_combine are used to match the sources
    in the image and the reference image, as well as make cuts based
    and other source properties. These values will likely have to
    adjusted based on your data.
    """
    img_dither = img.split('.')[1][0]+'_'
    ref_dither = ref.split('.')[1][0]+'_'

    img_sources = odi.sourcepath+img_dither+filter+'.allsource'
    ref_sources = odi.sourcepath+ref_dither+filter+'.allsource'

    x_img, y_img, id_img,ra_icrs_centroid_img,dec_icrs_centroid_img,source_sum_img,max_value_img,elongation_img, MAG_img, MERR_img, SKY_img, SERR_img, RAPERT_img, XPOS_img, YPOS_img, fwhm_img,peak_img = np.loadtxt(
        img_sources,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),unpack=True)

    x_ref, y_ref, id_ref,ra_icrs_centroid_ref,dec_icrs_centroid_ref,source_sum_ref,max_value_ref,elongation_ref, MAG_ref, MERR_ref, SKY_ref, SERR_ref, RAPERT_ref, XPOS_ref, YPOS_ref, fwhm_ref,peak_ref = np.loadtxt(
        ref_sources,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),unpack=True)

    img_catalog = SkyCoord(ra = ra_icrs_centroid_img*u.degree, dec = dec_icrs_centroid_img*u.degree)

    ref_catalog = SkyCoord(ra = ra_icrs_centroid_ref*u.degree, dec = dec_icrs_centroid_ref*u.degree)

    #id_img, d2d_img, d3d_img = ref_catalog.match_to_catalog_sky(img_catalog)

    #id_ref, d2d_ref, d3d_ref = img_catalog.match_to_catalog_sky(ref_catalog)

    id_img, id_ref, d2d, d3d = ref_catalog.search_around_sky(img_catalog,0.0001*u.deg)

    print len(id_img),len(id_ref)

    MAG_img    = MAG_img[id_img]
    MERR_img   = MERR_img[id_img]
    SKY_img    = SKY_img[id_img]
    SERR_img   = SERR_img[id_img]
    RAPERT_img = RAPERT_img[id_img]
    XPOS_img   = XPOS_img[id_img]
    YPOS_img   = YPOS_img[id_img]
    fwhm_img   = fwhm_img[id_img]
    peak_img   = peak_img[id_img]
    ra_icrs_centroid_img  =  ra_icrs_centroid_img[id_img]
    dec_icrs_centroid_img =  dec_icrs_centroid_img[id_img]

    MAG_ref    = MAG_ref[id_ref]
    MERR_ref   = MERR_ref[id_ref]
    SKY_ref    = SKY_ref[id_ref]
    SERR_ref   = SERR_ref[id_ref]
    RAPERT_ref = RAPERT_ref[id_ref]
    XPOS_ref   = XPOS_ref[id_ref]
    YPOS_ref   = YPOS_ref[id_ref]
    fwhm_ref   = fwhm_ref[id_ref]
    peak_ref   = peak_ref[id_ref]
    ra_icrs_centroid_ref  =  ra_icrs_centroid_ref[id_ref]
    dec_icrs_centroid_ref =  dec_icrs_centroid_ref[id_ref]


    #keep = np.where((SKY_img>0.0) & (SKY_ref > 0.0) & (peak_img>200) & (peak_ref >200.0) & (45000.0>peak_img) & (45000.0>peak_ref) & (peak_img < 100) & (peak_ref < 100))
    keep = np.where((np.array(peak_img)>1000.0) & (np.array(peak_ref) >1000.0)&(np.array(peak_img)<45000.0) & (np.array(peak_ref) <45000.0)
                    & (np.array(fwhm_img)<900.0) & (np.array(fwhm_ref) <900.0))

    magA =  np.array(MAG_img[keep[0]])
    magRef =  np.array(MAG_ref[keep[0]])
    raA, decA = ra_icrs_centroid_img[keep], dec_icrs_centroid_img[keep]
    raRef, decRef = ra_icrs_centroid_ref[keep], dec_icrs_centroid_ref[keep]

    with open('scale_stars.pos','w+') as f:
        for i,m in enumerate(magA):
            print >> f, raA[i], decA[i], raRef[i], decRef[i], magA[i], magRef[i]
            
    rat = np.power(10.0,-0.4*(magA-magRef))/1.0

    #print np.mean(rat),np.std(rat),len(rat)
    sigThreshold = 0.009
    n = 1
    sigTest = np.std(rat)
    if sigTest <= sigThreshold:
            scale = np.mean(rat)
            std = np.std(rat)
    else:
        while sigTest > sigThreshold:
            magTempA = magA
            magTempRef = magRef
            magA = magTempA[np.where(abs(rat-np.median(rat))<sigTest)]
            magRef = magTempRef[np.where(abs(rat-np.median(rat))<sigTest)]
            rat = np.power(10.0,-0.4*(magA-magRef))/1.0
            #for i,r in enumerate(rat):
            #print magA[i], magRef[i], r
            sigTest = np.std(rat)
            n = n + 1
            if n > 10:
                print "Iteration did not converge to sigma <", repr(sigThreshold),"for", img
                print "Quitting..."
                exit()
            #print len(rat), np.mean(rat), np.median(rat), np.std(rat), n
            #scale[img] = np.mean(rat)
            #std[img] = np.std(rat)
        scale = np.mean(rat)
        std = np.std(rat)
    return scale,std,len(rat)

def sdss_source_props_ota(img,ota):
    """
    Use photutils to get the elongation of all of the sdss sources
    can maybe use for point source filter
    Also fit a gaussian along a row and col of pixels passing
    through the center of the star
    """

    image = odi.reprojpath+'reproj_'+ota+'.'+str(img[16:])
    hdulist = odi.fits.open(image)
    data = hdulist[0].data

    sdss_source_file = odi.coordspath+'reproj_'+ota+'.'+str(img[16:-5])+'.sdssxy'

    x,y,ra,dec,g,g_err,r,r_err = np.loadtxt(sdss_source_file,usecols=(0,1,2,3,
                                                                      6,7,8,9),unpack=True)

    box_centers = zip(y,x)
    box_centers = np.reshape(box_centers,(len(box_centers),2))
    source_dict = {}
    total_fwhm = []
    for i,center in enumerate(box_centers):
        x1 = center[0]-50
        x2 = center[0]+50
        y1 = center[1]-50
        y2 = center[1]+50

        #print x1,x2,y1,y2,center
        box = data[x1:x2,y1:y2]
        col = data[x1:x2,int(center[1]-0.5):int(center[1]+0.5)]
        row = data[int(center[0]-0.5):int(center[0]+0.5),y1:y2]
        row = np.squeeze(row) - np.median(row)
        col = np.squeeze(col) - np.median(col)
        g_init = models.Gaussian1D(amplitude=250., mean=50, stddev=2.)
        fit_g = fitting.LevMarLSQFitter()
        pix = np.linspace(0,100,num=100)
        g_row = fit_g(g_init, pix, row)
        g_col = fit_g(g_init, pix, col)
        mean_fwhm = 0.5*(g_row.stddev*2.355+g_col.stddev*2.355)
        total_fwhm.append(mean_fwhm)
        #odi.plt.imshow(box)
        #odi.plt.plot(row)
        #odi.plt.plot(pix,g(pix))
        #plt.imshow(row2)
        #plt.show()
        mean, median, std = odi.sigma_clipped_stats(box, sigma=3.0)
        threshold = median + (std * 2.)
        segm_img = odi.detect_sources(box, threshold, npixels=20)
        source_props = odi.source_properties(box,segm_img)
        if len(source_props) > 0:
            columns = ['xcentroid', 'ycentroid','elongation','semimajor_axis_sigma','semiminor_axis_sigma']
            if i == 0:
                source_tbl = odi.properties_table(source_props,columns=columns)
            else:
                source_tbl.add_row((source_props[0].xcentroid,source_props[0].ycentroid,
                                    source_props[0].elongation,source_props[0].semimajor_axis_sigma,
                                    source_props[0].semiminor_axis_sigma))
    elong_med,elong_std = np.median(source_tbl['elongation']),np.std(source_tbl['elongation'])
    hdulist.close()
    return elong_med,elong_std,np.mean(total_fwhm),np.std(total_fwhm)
