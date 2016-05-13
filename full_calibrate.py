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

def full_sdssmatch(img1,img2,inst,gmaglim=19):
    """
    Fetch sdss images in each final stacked image
    Then match the sources in common between each
    image based on ra and dec.
    """
    odi.sdss_coords_full(img1,inst,gmaglim=gmaglim)
    img1_sdss_cat = img1[:-5]+'.sdssxy'
    img1_match = img1[:-5]+'.match.sdssxy'
    odi.sdss_coords_full(img2,inst,gmaglim=gmaglim)
    img2_sdss_cat = img2[:-5]+'.sdssxy'
    img2_match = img2[:-5]+'.match.sdssxy'
    
    x_1, y_1, ras_1,decs_1,psfMag_u_1,psfMagErr_u_1,psfMag_g_1,psfMagErr_g_1,psfMag_r_1,psfMagErr_r_1,psfMag_i_1,psfMagErr_i_1,psfMag_z_1,psfMagErr_z_1 = np.loadtxt(img1_sdss_cat,
                                                                                                                                                                     usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13),
                                                                                                                                                                     unpack=True)
    
    x_2, y_2, ras_2,decs_2,psfMag_u_2,psfMagErr_u_2,psfMag_g_2,psfMagErr_g_2,psfMag_r_2,psfMagErr_r_2,psfMag_i_2,psfMagErr_i_2,psfMag_z_2,psfMagErr_z_2 = np.loadtxt(img2_sdss_cat,
                                                                                                                                                                     usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13),
                                                                                                                                                                     unpack=True)
    

    img1_catalog = SkyCoord(ra = ras_1*u.degree, dec= decs_1*u.degree)
    
    img2_catalog = SkyCoord(ra = ras_2*u.degree, dec= decs_2*u.degree)
    
    id_img1, id_img2, d2d, d3d = img2_catalog.search_around_sky(img1_catalog,0.0001*u.deg)
    
    x_1              = x_1[id_img1]              
    y_1              = y_1[id_img1]              
    ras_1            = ras_1[id_img1]            
    decs_1           = decs_1[id_img1]           
    psfMag_u_1       = psfMag_u_1[id_img1]       
    psfMagErr_u_1    = psfMagErr_u_1[id_img1]    
    psfMag_g_1       = psfMag_g_1[id_img1]       
    psfMagErr_g_1    = psfMagErr_g_1[id_img1]    
    psfMag_r_1       = psfMag_r_1[id_img1]       
    psfMagErr_r_1    = psfMagErr_r_1[id_img1]    
    psfMag_i_1       = psfMag_i_1[id_img1]       
    psfMagErr_i_1    = psfMagErr_i_1[id_img1]    
    psfMag_z_1       = psfMag_z_1[id_img1]       
    psfMagErr_z_1    = psfMagErr_z_1[id_img1]
    
    img1_match_dict = OrderedDict([('x_1',x_1),('y_1',y_1),('ras_1',ras_1),
                                   ('decs_1',decs_1),('psfMag_u_1',psfMag_u_1),
                                   ('psfMagErr_u_1',psfMagErr_u_1),
                                   ('psfMag_g_1',psfMag_g_1),('psfMagErr_g_1',psfMagErr_g_1),
                                   ('psfMag_r_1',psfMag_r_1),('psfMagErr_r_1',psfMagErr_r_1),
                                   ('psfMag_i_1',psfMag_i_1),('psfMagErr_i_1',psfMagErr_i_1),
                                   ('psfMag_z_1',psfMag_z_1),('psfMagErr_z_1',psfMagErr_z_1)])

    img1_match_df = pd.DataFrame.from_dict(img1_match_dict)
    img1_match_df.to_csv(img1_match,index=False,sep= ' ',header=False)

    x_2              = x_2[id_img2]              
    y_2              = y_2[id_img2]              
    ras_2            = ras_2[id_img2]            
    decs_2           = decs_2[id_img2]           
    psfMag_u_2       = psfMag_u_2[id_img2]       
    psfMagErr_u_2    = psfMagErr_u_2[id_img2]    
    psfMag_g_2       = psfMag_g_2[id_img2]       
    psfMagErr_g_2    = psfMagErr_g_2[id_img2]    
    psfMag_r_2       = psfMag_r_2[id_img2]       
    psfMagErr_r_2    = psfMagErr_r_2[id_img2]    
    psfMag_i_2       = psfMag_i_2[id_img2]       
    psfMagErr_i_2    = psfMagErr_i_2[id_img2]    
    psfMag_z_2       = psfMag_z_2[id_img2]       
    psfMagErr_z_2    = psfMagErr_z_2[id_img2] 
    
    
    img2_match_dict = OrderedDict([('x_2',x_2),('y_2',y_2),('ras_2',ras_2),
                                   ('decs_2',decs_2),('psfMag_u_2',psfMag_u_2),
                                   ('psfMagErr_u_2',psfMagErr_u_2),
                                   ('psfMag_g_2',psfMag_g_2),('psfMagErr_g_2',psfMagErr_g_2),
                                   ('psfMag_r_2',psfMag_r_2),('psfMagErr_r_2',psfMagErr_r_2),
                                   ('psfMag_i_2',psfMag_i_2),('psfMagErr_i_2',psfMagErr_i_2),
                                   ('psfMag_z_2',psfMag_z_2),('psfMagErr_z_2',psfMagErr_z_2)])
    
    img2_match_df = pd.DataFrame.from_dict(img2_match_dict)
    img2_match_df.to_csv(img2_match,index=False,sep= ' ',header=False)
    
    return img1_match_df, img2_match_df

def sdss_source_props_full(img):
    """
    Use photutils to get the elongation of all of the sdss sources
    can maybe use for point source filter
    """
    hdulist = odi.fits.open(img)
    data = hdulist[0].data
    
    sdss_source_file = img[:-5]+'.match.sdssxy'
    
    x,y,ra,dec,g,g_err,r,r_err = np.loadtxt(sdss_source_file,usecols=(0,1,2,3,
                                                                      6,7,8,9),unpack=True)
    
    box_centers = zip(y,x)
    box_centers = np.reshape(box_centers,(len(box_centers),2))
    source_dict = {}
    for i,center in enumerate(box_centers):
        x1 = center[0]-50
        x2 = center[0]+50
        y1 = center[1]-50
        y2 = center[1]+50
        
        #print x1,x2,y1,y2,center
        box = data[x1:x2,y1:y2]
        #odi.plt.imshow(box)
        #plt.show()
        mean, median, std = odi.sigma_clipped_stats(box, sigma=3.0)
        threshold = median + (std * 2.)
        segm_img = odi.detect_sources(box, threshold, npixels=20)
        source_props = odi.source_properties(box,segm_img)
        columns = ['xcentroid', 'ycentroid','elongation','semimajor_axis_sigma','semiminor_axis_sigma']
        if i == 0:
            source_tbl = odi.properties_table(source_props,columns=columns)
        else:
            source_tbl.add_row((source_props[0].xcentroid,source_props[0].ycentroid, 
                                source_props[0].elongation,source_props[0].semimajor_axis_sigma,
                                source_props[0].semiminor_axis_sigma))
    elong_med,elong_std = np.median(source_tbl['elongation']),np.std(source_tbl['elongation'])
    hdulist.close()
    return elong_med,elong_std

def read_proc(file,filter):
    """
    Read the derived properties file to get useful values
    need for phot and dao find
    """
    filter_str = np.loadtxt(file,usecols=(2,),unpack=True,dtype=str)
    fwhm,bg_mean,bg_med,bg_std = np.loadtxt(file,usecols=(3,6,7,8),unpack=True)
    
    median_fwhm = np.median(fwhm[np.where(filter_str == filter)])
    median_bg_mean = np.median(bg_mean[np.where(filter_str == filter)])
    median_bg_median = np.median(bg_med[np.where(filter_str == filter)])
    median_bg_std = np.median(bg_std[np.where(filter_str == filter)])
    
    return median_fwhm,median_bg_mean,median_bg_median,median_bg_std

def get_airmass(image_list):
    """
    Calculate the median arimass of
    all the dithers in a given
    filter
    """
    airmasses = []
    for img in image_list:
        hdulist = odi.fits.open(img)
        airmasses.append(hdulist[0].header['airmass'])
        hdulist.close()
    return np.median(airmasses)

def calc_airmass():
    from pyraf import iraf
    if not os.path.isfile('setairmass.done'):
	iraf.astutil.setairmass.setParam('images', "msc*fits")          # Input images
	iraf.astutil.setairmass.setParam('intype', "beginning")    # Input keyword time stamp
	iraf.astutil.setairmass.setParam('outtype', "effective")    # Output airmass time stamp\n
	iraf.astutil.setairmass.setParam('ra', "ra")           # Right acsension keyword (hours)
	iraf.astutil.setairmass.setParam('dec', "dec")          # Declination keyword (degrees)
	iraf.astutil.setairmass.setParam('equinox', "radeceq")        # Equinox keyword (years)
	iraf.astutil.setairmass.setParam('st', "st")           # Local siderial time keyword (hours)
	iraf.astutil.setairmass.setParam('ut', "time-obs")     # Universal time keyword (hours)
	iraf.astutil.setairmass.setParam('date', "date-obs")     # Observation date keyword
	iraf.astutil.setairmass.setParam('exposure', "exptime")      # Exposure time keyword (seconds)
	iraf.astutil.setairmass.setParam('airmass', "airmass")      # Airmass keyword (output)
	iraf.astutil.setairmass.setParam('utmiddle', "utmiddle")     # Mid-observation UT keyword (output)
	iraf.astutil.setairmass.setParam('scale', 750.)           # The atmospheric scale height\n
	iraf.astutil.setairmass.setParam('show', 'yes')            # Print the airmasses and mid-UT?
	iraf.astutil.setairmass.setParam('update', 'yes')            # Update the image header?
	iraf.astutil.setairmass.setParam('override', 'yes')            # Override previous assignments?
	iraf.astutil.setairmass()
	with open('setairmass.done','w+') as f1:
	    print >> f1, True
    else:
	print 'setairmass already done'
    
def sdss_phot_full(img,fwhm,airmass):
    """
    Run phot on the matches sdss stars
    """
    from pyraf import iraf
    iraf.ptools(_doprint=0)
    # first grab the header and hang on to it so we can use other values
    hdulist = odi.fits.open(img)
    hdr1 = hdulist[0].header
    filter = hdr1['filter']
    hdulist.close()
    
    iraf.unlearn(iraf.phot,iraf.datapars,iraf.photpars,iraf.centerpars,iraf.fitskypars)
    iraf.apphot.phot.setParam('interactive',"no")
    iraf.apphot.phot.setParam('verify',"no")
    iraf.datapars.setParam('datamax',50000.)
    iraf.datapars.setParam('gain',"gain")
    iraf.datapars.setParam('ccdread','rdnoise')
    iraf.datapars.setParam('exposure',"exptime")

    iraf.datapars.setParam('filter',"filter")
    iraf.datapars.setParam('obstime',"time-obs")
    iraf.datapars.setParam('sigma',"INDEF")
    iraf.photpars.setParam('zmag',0.)
    iraf.centerpars.setParam('cbox',9.)
    iraf.centerpars.setParam('maxshift',3.)
    iraf.fitskypars.setParam('salgorithm',"median")
    iraf.fitskypars.setParam('dannulus',10.)
    
    if not os.path.isfile(img[0:-5]+'.sdssphot'): # only do this once
        print 'phot-ing the sdss sources in ', filter
        iraf.datapars.setParam('xairmass',float(airmass))
        iraf.datapars.setParam('fwhmpsf',float(fwhm))
        iraf.photpars.setParam('apertures',5.*float(fwhm)) # use a big aperture for this
        iraf.fitskypars.setParam('annulus',6.*float(fwhm))
        iraf.apphot.phot(image=img, coords=img[:-5]+'.match.sdssxy', output=img[0:-5]+'.phot.1')
        phot_tbl = img[0:-5]+'.sdssphot'
        with open(phot_tbl,'w+') as txdump_out :
            iraf.ptools.txdump(textfiles=img[0:-5]+'.phot.1', fields="id,mag,merr,msky,stdev,rapert,xcen,ycen,ifilter,xairmass,image",expr='yes', headers='no', Stdout=txdump_out)
    
        outputfile_clean = open(phot_tbl.replace('.sdssphot','_clean.sdssphot'),"w")
        for line in open(phot_tbl,"r"):
            if not 'INDEF' in line:
                outputfile_clean.write(line)
            if 'INDEF' in line:
                outputfile_clean.write(line.replace('INDEF','999'))
        outputfile_clean.close()
        os.rename(phot_tbl.replace('.sdssphot','_clean.sdssphot'),phot_tbl)


def getfwhm_full_sdss(img, radius=4.0, buff=7.0, width=5.0):
    '''
    Get a fwhm estimate for the image using the SDSS catalog stars and IRAF imexam (SLOW, but works)
    Adapted from Kathy's getfwhm script (this implementation is simpler in practice)
    '''
    coords = img[:-5]+'.match.sdssxy'
    outputfile = img[0:-5]+'.sdssmatch.fwhm.log'

    iraf.tv.rimexam.setParam('radius',radius)
    iraf.tv.rimexam.setParam('buffer',buff)
    iraf.tv.rimexam.setParam('width',width)
    iraf.tv.rimexam.setParam('rplot',20.)
    iraf.tv.rimexam.setParam('center','yes')
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
    peak,gfwhm = np.loadtxt(outputfile, usecols=(9,10), unpack=True)
    
    return peak,gfwhm


        
def apcor_sdss(img,fwhm,inspect=False):
    """
    Use csv tables produced by full_sdssmatch
    """
    from pyraf import iraf
    import matplotlib.pyplot as plt
    from scipy import interpolate
    from matplotlib.colors import LogNorm
    from astropy.visualization import *
    from astropy.visualization.mpl_normalize import ImageNormalize
    iraf.ptools(_doprint=0)
    sdss_source_file = img[:-5]+'.match.sdssxy'
    
    sdss_phot_file = img[0:-5]+'.sdssphot'
    
    sdss_MAG, sdss_MERR, sdss_SKY, sdss_SERR, sdss_RAPERT, sdss_XPOS, sdss_YPOS = np.loadtxt(img[0:-5]+'.sdssphot', 
                                                                                             usecols=(1,2,3,4,5,6,7), 
                                                                                             dtype=float, unpack=True)
    
    x,y,ra,dec,g,g_err,r,r_err = np.loadtxt(sdss_source_file,usecols=(0,1,2,3,
                                                                      6,7,8,9),unpack=True)
    peak,gfwhm = getfwhm_full_sdss(img)
    
    #plt.clf()
    #plt.hist(peak[np.where(peak > 12000.0)])
    #plt.show()
    
    aps = []
    for i in np.arange(1,7.5,0.5):
        aps.append(fwhm*i)
    aps_str = str('"'+repr(aps[0])+','+repr(aps[1])+','
                  +repr(aps[2])+','+repr(aps[3])+','
                  +repr(aps[4])+','+repr(aps[5])+','
                  +repr(aps[6])+','+repr(aps[7])+','
                  +repr(aps[8])+','+repr(aps[9])+','
                  +repr(aps[10])+','+repr(aps[11])+','
                  +repr(aps[12])+'"')

    hdulist = odi.fits.open(img)
    hdr1 = hdulist[0].header
    data = hdulist[0].data
    filter = hdr1['filter']
    hdulist.close()
    
    iraf.unlearn(iraf.phot,iraf.datapars,iraf.photpars,iraf.centerpars,iraf.fitskypars)
    iraf.apphot.phot.setParam('interactive',"no")
    iraf.apphot.phot.setParam('verify',"no")
    iraf.datapars.setParam('datamax',50000.)
    iraf.datapars.setParam('gain',"gain")
    iraf.datapars.setParam('ccdread','rdnoise')
    iraf.datapars.setParam('exposure',"exptime")
    #iraf.datapars.setParam('itime',2700.0)

    iraf.datapars.setParam('filter',"filter")
    iraf.datapars.setParam('obstime',"time-obs")
    iraf.datapars.setParam('sigma',"INDEF")
    iraf.photpars.setParam('zmag',0.)
    iraf.centerpars.setParam('cbox',9.)
    iraf.centerpars.setParam('maxshift',3.)
    iraf.fitskypars.setParam('salgorithm',"median")
    iraf.fitskypars.setParam('dannulus',10.)
    
    phot_tbl = img[0:-5]+'.apcor'
    if not os.path.isfile(phot_tbl):
        print 'running phot over', aps_str
        iraf.datapars.setParam('fwhmpsf',fwhm)
        iraf.photpars.setParam('apertures',aps_str)
        iraf.fitskypars.setParam('annulus',6.5*fwhm)
        iraf.apphot.phot(image=img, coords=sdss_source_file, output=img[0:-5]+'.apcor.1')
        with open(phot_tbl,'w+') as txdump_out :
            iraf.ptools.txdump(textfiles=img[0:-5]+'.apcor.1', fields="ID,RAPERT,XCEN,YCEN,FLUX,MAG,MERR", expr="yes", headers='no', Stdout=txdump_out)
        txdump_out.close()
        outputfile_clean = open(phot_tbl.replace('.apcor','_clean.apcor'),"w")
        for line in open(phot_tbl,"r"):
            if not 'INDEF' in line:
                outputfile_clean.write(line)
            if 'INDEF' in line:
                outputfile_clean.write(line.replace('INDEF','999'))
        outputfile_clean.close()
        os.rename(phot_tbl.replace('.apcor','_clean.apcor'),phot_tbl)
        
    peak_top1per = 49000.0
    star_flux = {}
    star_mag_diff = {}
    star_mag_diff1x = {}
    star_mag = {}
    star_positions = {}
    for i,line in enumerate(open(img[0:-5]+'.apcor',"r")):
        flux = [float(x) for x in line.split()[16:29]]
        mag = [float(x) for x in line.split()[29:42]]
        err = [float(x) for x in line.split()[42:55]]
        position = [float(x) for x in line.split()[14:16]]
        if (peak[i] >= 18000.0 and peak[i] <= 45000.0  and 
            (np.abs(gfwhm[i] - np.median(gfwhm[np.where(gfwhm < 20.0)])) < np.std(gfwhm[np.where(gfwhm < 20.0)]))
            and np.max(mag) != 999.0 and g[i] <= 17):
            if inspect == True:
                center = position
                x1 = center[0]-75
                x2 = center[0]+75
                y1 = center[1]-75
                y2 = center[1]+75
                box = data[y1:y2,x1:x2]
                plt.figure()
                norm = ImageNormalize(stretch=LogStretch())
                odi.plt.imshow(box,norm=norm)
                max_count = np.max(box)
                odi.plt.title('max counts ='+str(max_count))
                plt.show()
                star_check = raw_input('Use star for ap correction: (y/n) ')
                if star_check == 'y':
                    star_flux[i] = flux
                    star_mag[i] = mag
                    star_positions[i] = position
            else:
                star_positions[i] = position
                star_flux[i] = flux
                star_mag[i] = mag
    for key in star_mag:
        diffs = []
        diffs1x = []
        for m in range(len(star_mag[key])-1):
            diffs.append(star_mag[key][m+1] - star_mag[key][m])
            diffs1x.append(star_mag[key][m] - star_mag[key][0])
        star_mag_diff[key] = diffs
        star_mag_diff1x[key] = diffs1x
    combine_mag_diffs = []
    combine_mag_diffs1x = []
    for key in star_mag_diff:
        x = np.arange(1,7.0,0.5)
        y = star_mag_diff[key]
        z = star_mag_diff1x[key]
        #plt.plot(x,y,'o')
        #plt.show()
        combine_mag_diffs.append(y)
        combine_mag_diffs1x.append(z)
        #tck = interpolate.splrep(x, y, s=0)
        #xnew = np.arange(1,6,0.25)
        #ynew = interpolate.splev(xnew,tck,der=1)
        #print len(np.arange(1,6.5,0.5)),len(star_mag_diff[key])
        #plt.plot(np.arange(1,6.5,0.5),star_mag_diff[key],'o')
        #plt.plot(xnew,ynew,'-')
        #plt.axhline(y=0)
        #plt.show()
    combine_mag_diffs = np.reshape(combine_mag_diffs,(len(star_mag_diff.keys()),len(x)))
    combine_mag_diffs_med = []
    combine_mag_diffs_std = []
    combine_mag_diffs1x = np.reshape(combine_mag_diffs1x,(len(star_mag_diff1x.keys()),len(x)))
    combine_mag_diffs1x_mean = []
    combine_mag_diffs1x_med = []
    combine_mag_diffs1x_std = []
    combine_mag_diffs1x_sem = []
    for j in range(len(x)):
        if j !=0 :
            star_test = combine_mag_diffs1x[:,j]
            sig_test = np.std(combine_mag_diffs1x[:,j])
            #print sig_test,len(star_test)
            med_test = np.median(combine_mag_diffs1x[:,j])
            keep = star_test[np.where( np.abs(star_test-med_test) < 0.010)]
            #plt.hist(keep)
            #plt.show()
            combine_mag_diffs1x_mean.append(np.mean(keep))
            combine_mag_diffs1x_med.append(np.median(keep))
            combine_mag_diffs1x_std.append(np.std(keep))
            combine_mag_diffs1x_sem.append(np.std(keep)/np.sqrt(len(keep)))
            if j == 7:
                print 'using ',len(keep), 'stars in ap correction'
        else:
            combine_mag_diffs1x_mean.append(np.mean(combine_mag_diffs1x[:,j]))
            combine_mag_diffs1x_med.append(np.median(combine_mag_diffs1x[:,j]))
            combine_mag_diffs1x_std.append(np.std(combine_mag_diffs1x[:,j]))
            combine_mag_diffs1x_sem.append(np.std(combine_mag_diffs1x[:,j])/np.sqrt(len(combine_mag_diffs1x[:,j])))
        
        combine_mag_diffs_med.append(np.median(combine_mag_diffs[:,j]))
        combine_mag_diffs_std.append(np.std(combine_mag_diffs[:,j]))
    
    tck = interpolate.splrep(x, combine_mag_diffs_med, s=0)
    xnew = np.arange(1,7.0,0.25)
    ynew = interpolate.splev(xnew,tck,der=0)
    ynew_der = interpolate.splev(xnew,tck,der=1)
    plt.figure(figsize=(14,7))
    ax1 = plt.subplot(121)
    ax2 = plt.subplot(122)
    ax1.errorbar(x,combine_mag_diffs_med,yerr=combine_mag_diffs_std,fmt='o',label='Median mag. Diff.')
    ax1.plot(xnew,ynew,'r-',label='Spline fit')
    ax1.set_ylabel('Mag(n+1) - Mag(n)')
    ax1.set_xlabel('n $ \cdot $ fwhm')
    ax2.plot(xnew,ynew_der,'r-')
    ax2.axhline(y=0)
    ax2.set_ylabel('Value of Spline 1st Derivative')
    ax2.set_xlabel('n $ \cdot $ fwhm')
    ax1.legend(loc=4)
    plt.tight_layout()
    plt.show()
    
    # print combine_mag_diffs1x_mean
    # print combine_mag_diffs1x_std
    # print combine_mag_diffs1x_sem
    
    for i in range(len(x)):
        plt.errorbar(x[i],combine_mag_diffs1x_mean[i],yerr=combine_mag_diffs1x_std[i],fmt='bo')
    plt.ylim(-0.45,0.1)
    plt.xlim(-0.1,7.5)
    plt.xlabel('n fwhm')
    plt.ylabel('ap correction')
    plt.show()
    
    apcor = np.median(combine_mag_diffs1x_mean[5:10:1])
    apcor_std = combine_mag_diffs1x_std[7]
    apcor_sem = combine_mag_diffs1x_sem[7]
    apcor_med = combine_mag_diffs1x_med[7]
    
    print 'aperture corr. = {0:6.3f}'.format(apcor)
    print 'aperture corr. med. = {0:6.3f}'.format(apcor_med)
    print 'aperture corr. std = {0:6.3f}'.format(apcor_std)
    print 'aperture corr. sem = {0:6.3f}'.format(apcor_sem)
    
    return apcor, apcor_std, apcor_sem
        
def calibrate_match(img1, img2, fwhm1, fwhm2, airmass1, airmass2):
    """
    Solve color equations
    """
    try:
        from pyraf import iraf
        from astropy.io import fits
        import numpy as np
        from scipy import stats
        import scipy.optimize as opt
        import matplotlib.pyplot as plt
    except ImportError:
        print 'You need some non-core python packages and a working IRAF to run this program'
        print "Try 'pip install astropy numpy scipy matplotlib pyraf' and try again"

    img_root = img1[:-7]
    
    # values determined by ralf/daniel @ wiyn
    kg = 0.20
    kr = 0.12
    ki = 0.058

    # you're going to need the average stellar fwhm to compute a aperture size
    # ralf or steven probably write one to the image header during QR/etc
    # just use that value here

    # first grab the header and hang on to it so we can use other values
    hdulist = fits.open(img1)
    hdr1 = hdulist[0].header
    hdulist.close()

    # for both images
    hdulist = fits.open(img2)
    hdr2 = hdulist[0].header
    hdulist.close()

    # read in the phot output as a string because we need to get rid of the indefs
    gMAG, gMERR, gSKY, gSERR, gRAPERT, gXPOS, gYPOS = np.loadtxt(img1[0:-5]+'.sdssphot', usecols=(1,2,3,4,5,6,7), dtype=float, unpack=True)
    iMAG, iMERR, iSKY, iSERR, iRAPERT, iXPOS, iYPOS = np.loadtxt(img2[0:-5]+'.sdssphot', usecols=(1,2,3,4,5,6,7), dtype=float, unpack=True)

    gXAIRMASS = np.loadtxt(img1[0:-5]+'.sdssphot', usecols=(9,), dtype=str, unpack=True)
    iXAIRMASS = np.loadtxt(img2[0:-5]+'.sdssphot', usecols=(9,), dtype=str, unpack=True)
    
    gFILTER = np.loadtxt(img1[0:-5]+'.sdssphot', usecols=(8,), dtype=str, unpack=True)
    iFILTER = np.loadtxt(img2[0:-5]+'.sdssphot', usecols=(8,), dtype=str, unpack=True)

    gID = np.loadtxt(img1[0:-5]+'.sdssphot', usecols=(0,), dtype=int, unpack=True)
    iID = np.loadtxt(img2[0:-5]+'.sdssphot', usecols=(0,), dtype=int, unpack=True)

    # keep the actual ID number to select from SDSS stars
    gID_keep = gID - 1
    iID_keep = iID - 1
    keep = list(set(gID_keep).intersection(iID_keep))

    # and keep the common elements between g and i using their list index
    keepg = [i for i,element in enumerate(gID) if element in iID]
    keepi = [i for i,element in enumerate(iID) if element in gID]

    # check to see if we're actually getting the same star across all of the files
    # for i in range(len(keep)):
    #     print keep[i]+1, gID[keepg[i]], iID[keepi[i]]
    # and how many
    # print len(keep), len(keepg), len(keepi)

    # read in the the SDSS catalog values
    x, y, ra, dec, u, ue, g, ge, r, re, i, ie, z, ze = np.loadtxt(img1[:-5]+'.match.sdssxy', usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13), unpack=True)

    # pick out the ones that match the good phot stars
    g, ge, r, re, i, ie = np.array(g[keep]), np.array(ge[keep]), np.array(r[keep]), np.array(re[keep]), np.array(i[keep]), np.array(ie[keep])

    # and reduce the other vectors
    gXPOS, gYPOS, gMAG, gMERR, gSKY, gSERR, iMAG, iMERR, iSKY, iSERR = np.array(gXPOS[keepg]), np.array(gYPOS[keepg]), np.array(gMAG[keepg]), np.array(gMERR[keepg]), np.array(gSKY[keepg]), np.array(gSERR[keepg]), np.array(iMAG[keepi]), np.array(iMERR[keepi]), np.array(iSKY[keepi]), np.array(iSERR[keepi])

    # keep the airmasses and aperture radii as single values
    if gXAIRMASS[0] != 'INDEF':
        gXAIRMASS, iXAIRMASS = gXAIRMASS.astype(float)[0], iXAIRMASS.astype(float)[0]
    else:
        gXAIRMASS, iXAIRMASS = airmass1, airmass2
    gRAPERT, iRAPERT = gRAPERT[0], iRAPERT[0]

    # apply airmass extinction correction to instrumental magnitudes
    g0 = gMAG - kg*gXAIRMASS
    if iFILTER[0].endswith('i'):
        print 'you gave me an i-band image, proceeding...'
        i0 = iMAG - ki*iXAIRMASS
        filterName = 'i'
        # determine catalog color and error
        gi = g - i
        gie = np.sqrt(ge**2 + ie**2)
    elif iFILTER[0].endswith('r'):
        print 'you gave me an r-band image, proceeding...'
        i0 = iMAG - kr*iXAIRMASS    
        filterName = 'r'
        # determine catalog color and error
        i = r
        ie = re
        gi = g - r
        gie = np.sqrt(ge**2 + re**2)

    # from here on, all i variables represent either i or r depending on what the user input
    # determine instrumental color and its associated error
    gi0 = g0 - i0
    giMERR = np.sqrt(gMERR**2 + iMERR**2)

    # find the difference between instrumental i or r and catalog value & error
    di = i - i0
    die = np.sqrt(ie**2 + iMERR**2)

    podicut, sdsscut = 0.01, 0.03
    print np.median(gSERR), np.median(iSERR)
    # cuts for better fits go here
    errcut = [j for j in range(len(gMERR)) if (gMERR[j] < podicut and iMERR[j] < podicut and ge[j] < sdsscut and ie[j] < sdsscut and gSKY[j] > np.median(gSERR) and iSKY[j] > np.median(iSERR))]

    with open('photcal_stars.pos','w+') as f1:
        for i, xp in enumerate(gXPOS[errcut]):
            print >> f1, xp, gYPOS[i]
            
    print len(gi0[errcut])

    # fit color term
    # linear lsq with numpy.polyfit
    p, pcov = np.polyfit(gi0[errcut], gi[errcut], 1, cov=True)
    perr = np.sqrt(np.diag(pcov))
    mu_gi, zp_gi, std_mu_gi, std_zp_gi = p[0], p[1], perr[0], perr[1]

    # print mu_gi, zp_gi, std_mu_gi, std_zp_gi

    # do a sigma clip based on the rms of the data from the first fit
    xplt1 = gi0[errcut]
    yplt1 = mu_gi*xplt1 + zp_gi

    dy1 = yplt1 - gi[errcut]

    # print std_zp_i
    # this actually pulls out the clipped values
    gi0_2 = np.array([col for j,col in enumerate(gi0[errcut]) if (abs(dy1[j]) < dy1.std())])
    gi_2 = np.array([col for j,col in enumerate(gi[errcut]) if (abs(dy1[j]) < dy1.std())])

    # linear lsq with numpy.polyfit
    p, pcov = np.polyfit(gi0_2, gi_2, 1, cov=True)
    perr = np.sqrt(np.diag(pcov))
    mu_gi, zp_gi, std_mu_gi, std_zp_gi = p[0], p[1], perr[0], perr[1]
    
    # set up 95% confidence interval calculation
    conf = 0.95
    alpha=1.-conf	# significance
    n=gi0_2.size	# data sample size
    x = np.arange(-1.0,3.5,0.025)
    # Auxiliary definitions
    mse=1./(n-2.)* np.sum((gi_2-(mu_gi*gi0_2 + zp_gi))**2)	# Scatter of data about the model (mean square error)
    stdev = np.sqrt(mse)
    sxd=np.sum((gi0_2-gi0_2.mean())**2) # standard deviation of data
    sx=(x-gi0_2.mean())**2	# fit residuals
    
    # Quantile of Student's t distribution for p=1-alpha/2
    q=stats.t.ppf(1.-alpha/2.,n-2)
    
    # 95% Confidence band
    dy=q*np.sqrt(mse*(1./n + sx/sxd ))
    mu_ucb=mu_gi*x + zp_gi +dy	# Upper confidence band
    mu_lcb=mu_gi*x + zp_gi -dy	# Lower confidence band


    print '--------------------------------------------------------------------------'
    print 'Here are the fit values:'
    print 'mu_g'+filterName+'      std_mu_g'+filterName+'  zp_g'+filterName+'      std_zp_g'+filterName
    print '{0:10.7f} {1:10.7f} {2:10.7f} {3:10.7f}'.format(mu_gi, std_mu_gi, zp_gi, std_zp_gi)

    # fit zero point
    # linear lsq with numpy.polyfit
    p, pcov = np.polyfit(gi[errcut], di[errcut], 1, cov=True)
    perr = np.sqrt(np.diag(pcov))
    eps_gi, zp_i, std_eps_gi, std_zp_i = p[0], p[1], perr[0], perr[1]

    # print eps_gi, zp_i, std_eps_gi, std_zp_i

    # do a sigma clip based on the rms of the data from the first fit
    xplt2 = gi[errcut]
    yplt2 = eps_gi*xplt2 + zp_i

    dy2 = yplt2 - di[errcut]

    # print std_zp_i
    # this actually pulls out the clipped values
    gi_3 = np.array([col for j,col in enumerate(gi[errcut]) if (abs(dy2[j]) < dy2.std())])
    di_3 = np.array([col for j,col in enumerate(di[errcut]) if (abs(dy2[j]) < dy2.std())])
    gX_3 = np.array([col for j,col in enumerate(gXPOS[errcut]) if (abs(dy2[j]) < dy2.std())])
    gY_3 = np.array([col for j,col in enumerate(gYPOS[errcut]) if (abs(dy2[j]) < dy2.std())])

    # linear lsq with numpy.polyfit
    p, pcov = np.polyfit(gi_3, di_3, 1, cov=True)
    perr = np.sqrt(np.diag(pcov))
    eps_gi, zp_i, std_eps_gi, std_zp_i = p[0], p[1], perr[0], perr[1]
    print 'eps_g'+filterName+'     std_eps_g'+filterName+' zp_'+filterName+'        std_zp_'+filterName
    print '{0:10.7f} {1:10.7f} {2:10.7f} {3:10.7f}'.format(eps_gi, std_eps_gi, zp_i, std_zp_i)
    
    #zp_check=[]
    #for i in [2,3,4]:
        #for j in [2,3,4]:
            #try:
                #zp_chk = ota_zp(gX_3, gY_3, gi_3, di_3, i, j)
                #zp_check.append(zp_chk)
            #except:
                #zp_check.append(0.0)
    #print np.std(np.array(zp_check))
    #print zp_check
    
    # set up 95% confidence interval calculation
    conf = 0.95
    alpha=1.-conf	# significance
    n=gi_3.size	# data sample size
    x = np.arange(-1.0,3.5,0.025)
    # Auxiliary definitions
    mse=1./(n-2.)* np.sum((di_3-(eps_gi*gi_3 + zp_i))**2)	# Scatter of data about the model (mean square error)
    stdev = np.sqrt(mse)
    sxd=np.sum((gi_3-gi_3.mean())**2) # standard deviation of data
    sx=(x-gi_3.mean())**2	# fit residuals
    
    # Quantile of Student's t distribution for p=1-alpha/2
    q=stats.t.ppf(1.-alpha/2.,n-2)
    
    # 95% Confidence band
    dy=q*np.sqrt(mse*(1./n + sx/sxd ))
    eps_ucb=eps_gi*x + zp_i +dy	# Upper confidence band
    eps_lcb=eps_gi*x + zp_i -dy	# Lower confidence band

    # make a diagnostic plot
    xplt = np.arange(-2,6,0.1)
    yplt = mu_gi*xplt + zp_gi

    plt.subplot(211)
    plt.scatter(gi0[errcut], gi[errcut], facecolor='red', edgecolor='none', s=3)
    plt.scatter(gi0_2, gi_2, facecolor='black', edgecolor='none', s=3)
    plt.plot(xplt, yplt, 'r-', lw=1, alpha=1, label='fit')
    # put 2xRMS on the plot
    plt.fill_between(x, mu_ucb, mu_lcb, facecolor='blue', edgecolor='none', alpha=0.2, label='2x RMS sigma clipping region')
    plt.xlim(-1,3.5)
    plt.xlabel('$g_0 - '+filterName+'_0$ (ODI)')
    plt.ylim(-1,3.5)
    plt.ylabel('$g - '+filterName+'$ (SDSS)')
    plt.text(-0.9, 3.0, '$\mu_{g'+filterName+'} = %.7f \pm %.7f$'%(mu_gi,std_mu_gi))
    plt.text(-0.9, 2.5, '$\mathrm{zp}_{g'+filterName+'} = %.7f \pm %.7f$'%(zp_gi,std_zp_gi))
    # plt.legend(loc=3)

    plt.subplot(212)
    xplt = np.arange(-2,6,0.1)
    yplt = eps_gi*xplt + zp_i
    # plt.plot([-2,-2],[0,0], 'k--')
    plt.scatter(gi[errcut], di[errcut], facecolor='red', edgecolor='none', s=3)
    plt.scatter(gi_3, di_3, facecolor='black', edgecolor='none', s=3)
    plt.plot(xplt, yplt, 'r-', lw=1, alpha=1, label='fit')
    plt.fill_between(x, eps_ucb, eps_lcb, facecolor='blue', edgecolor='none', alpha=0.2, label='2x RMS sigma clipping region')
    plt.xlim(-1,3.5)
    plt.ylim(zp_i+1.0,zp_i-1.0)
    plt.xlabel('$g - '+filterName+'$ (SDSS)')
    plt.ylabel('$'+filterName+' - '+filterName+'_0$ (SDSS - ODI)')
    plt.text(-0.9, zp_i-0.8, '$\epsilon_{g'+filterName+'} = %.5f \pm %.5f$'%(eps_gi,std_eps_gi))
    plt.text(-0.9, zp_i-0.6, '$\mathrm{zp}_{'+filterName+'} = %.5f \pm %.5f$'%(zp_i,std_zp_i))
    plt.tight_layout()
    plt.savefig(img_root+'_photcal.pdf')
    
    plt.clf()
    plt.scatter(gXPOS, gYPOS, c='red', edgecolor='none')
    plt.xlabel('X pixel')
    plt.ylabel('Y pixel')
    plt.xlim(0,13500)
    plt.ylim(0,13500)
    plt.savefig(img_root+'_photmap.pdf')
    
    # make a cmd of the ODI photometry of all the SDSS stars for reference
    g0 = gMAG - (kg*gXAIRMASS)
    i0 = iMAG - (ki*iXAIRMASS)
    gmi = mu_gi*(g0-i0) + zp_gi
    
    i_mag = i0 + eps_gi*gmi + zp_i #- cal_A_i 
    g_mag = gmi + i_mag
    
    plt.clf()
    plt.scatter(gmi, i_mag, c='red', s=3, edgecolor='none')
    plt.xlabel('$g-r$')
    plt.ylabel('$r$')
    plt.xlim(-1,2)
    plt.ylim(24,14)
    plt.savefig(img_root+'_photcmd.pdf')

    # print out a steven style help file, no writing to headers YET
    with open(img_root+'_help.txt','w+') as f1:
        print >> f1, "this has some information about the calibration. don't panic."
        print >> f1, ""
        print >> f1, "this is the revised (Feb 2015) version of pODI - SDSS calibrations"
        print >> f1, "   it is run on matched pairs of images (g+i, for UCHVC project)"
        print >> f1, ""
        print >> f1, "it follows the extremely standard method of photometric calibrations:"
        print >> f1, ""
        print >> f1, "g-i = mu_gi ( g0 - i0 ) + ZP_gi"
        print >> f1, "i = i0 + eps_gi ( g - i ) + ZP_i"
        print >> f1, ""
        print >> f1, "   where g0 = g_i - k_g * X_g  include airmass extinction"
        print >> f1, "         i0 = i_i - k_i * X_i"
        print >> f1, "Fits generate errors on mu/eps/ZP and also rms for both"
        print >> f1, ""
        print >> f1, "g_i/i_i are instrumental magnitudes, measured in apertures 5x FWHM"
        print >> f1, ""
        print >> f1, "all of these coefficients are saved to both g&i image headers,"
        print >> f1, "    and are reproduced below."
        print >> f1, ""
        print >> f1, "in particular, this is the calibration for $!gal"
        print >> f1, ""
        print >> f1, "  name          symbol     IMHEAD     value"
        print >> f1, "----------------------------------------------------"
        print >> f1, "  extn coeff      k_g      F_KG       {0:.7f}".format(kg)
        print >> f1, "  extn coeff      k_i      F_KI       {0:.7f}".format(ki)
        print >> f1, "  airmass in g    X_g      F_XG       {0:.7f}".format(gXAIRMASS)
        print >> f1, "  airmass in i    X_i      F_XI       {0:.7f}".format(iXAIRMASS)
        print >> f1, " - - - - - - - - - - - - - - - - - - - - - - - - - -"
        print >> f1, "  g-i color term  mu_gi    F_MU_GI    {0:.7f}".format(mu_gi)
        print >> f1, "  g-i c.t. err    mue_gi   F_MUE_GI   {0:.7f}".format(std_mu_gi)
        print >> f1, "  g-i zeropoint   ZP_gi    F_ZP_GI    {0:.7f}".format(zp_gi)
        print >> f1, "  g-i ZP err      ZPE_gi   F_ZPE_GI   {0:.7f}".format(std_zp_gi)
        print >> f1, "  g-i fit RMS     rms      F_RMS_GI   {0:.7f}".format(dy1.std())
        print >> f1, " - - - - - - - - - - - - - - - - - - - - - - - - - -"
        print >> f1, "  i color term    eps_gi   F_EPS_GI   {0:.7f}".format(eps_gi)
        print >> f1, "  i c.t. err      epse_gi  F_EPSE_GI  {0:.7f}".format(std_eps_gi)
        print >> f1, "  i zeropoint     ZP_i     F_ZP_I     {0:.7f}".format(zp_i)
        print >> f1, "  i ZP err        ZPe_i    F_ZPE_I    {0:.7f}".format(std_zp_i)
        print >> f1, "  i fit RMS       rms      F_RMS_I    {0:.7f}".format(dy2.std())
        print >> f1, "----------------------------------------------------"
        print >> f1, "other details:"
        print >> f1, "  FWHM PSF [px]   fwhm    FWHMPSF    [see header]"
        print >> f1, "  FWHM [arcsec] g fwhm    F_AVGSEE   {0:.5f}".format(0.11*gRAPERT/5)
        print >> f1, "  FWHM [arcsec] i fwhm    F_AVGSEE   {0:.5f}".format(0.11*iRAPERT/5)
        print >> f1, "  phot aperture (5xFWHM) g [arcsec]  {0:.5f}".format(0.11*gRAPERT)
        print >> f1, "  phot aperture (5xFWHM) i [arcsec]  {0:.5f}".format(0.11*iRAPERT)
        print >> f1, "----------------------------------------------------"
        print >> f1, "photometric error cuts:"
        print >> f1, "  maximum acceptable pODI PHOT error: {0:.4f}".format(podicut)
        print >> f1, "  maximum acceptable sdss phot error: {0:.4f}".format(sdsscut)
        print >> f1, "  N_stars surviving error cuts: {0:4d}".format(len(gi[errcut]))
        print >> f1, "    N_stars surviving sigma clip (i-i0 vs g-i plot): {0:4d}".format(len(gi_3))
    print '--------------------------------------------------------------------------'
    print 'Done! I saved some important information in the following files for you:'
    print 'SDSS raw catalog values (csv):         ', img_root+'.sdss'
    print 'SDSS catalog values w/ x,y positions:  ', img_root+'.sdssxy'
    print 'Instrumental ODI magnitudes per image: ', img_root+'*.sdssphot'
    print 'Calibration fit diagnostic plots:      ', img_root+'_photcal.pdf'
    print 'Final calibration values:              ', img_root+'_help.txt'


# How to use in odi_process
#g_img = 'GCPair-F1_odi_g.fits'
#r_img = 'GCPair-F1_odi_r.fits'

#odi.full_sdssmatch(g_img,r_img,inst)

#median_fwhm,median_bg_mean,median_bg_median,median_bg_std = odi.read_proc('derived_props.txt','odi_g')
#print median_fwhm,median_bg_mean,median_bg_median,median_bg_std

#airmass_g = odi.get_airmass(images_g)

#odi.sdss_phot(g_img,median_fwhm,airmass_g)

#median_fwhmr,median_bg_meanr,median_bg_medianr,median_bg_stdr = odi.read_proc('derived_props.txt','odi_r')
#print median_fwhmr,median_bg_meanr,median_bg_medianr,median_bg_stdr

#airmass_r = odi.get_airmass(images_r)

#odi.sdss_phot(r_img,median_fwhmr,airmass_r)

#odi.calibrate_match(g_img,r_img,median_fwhm,median_fwhmr,airmass_g,airmass_r)