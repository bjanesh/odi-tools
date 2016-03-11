#! /usr/local/bin/python
def escut(image, pos_file, fwhm, peak):
    # input image file name, file name with matched source positions, **np.array of fwhm measurements for each source
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import stats
    from pyraf import iraf
    # import sewpy
    import os
    from matplotlib.path import Path
    
    iraf.images(_doprint=0)
    iraf.tv(_doprint=0)
    iraf.ptools(_doprint=0)
    iraf.noao(_doprint=0)
    iraf.digiphot(_doprint=0)
    iraf.photcal(_doprint=0)
    iraf.apphot(_doprint=0)  
    iraf.imutil(_doprint=0)
    
    iraf.unlearn(iraf.phot,iraf.datapars,iraf.photpars,iraf.centerpars,iraf.fitskypars)
    iraf.apphot.phot.setParam('interactive',"no")
    iraf.apphot.phot.setParam('verify',"no")
    iraf.datapars.setParam('datamax',50000.)
    iraf.datapars.setParam('gain',"gain")
    iraf.datapars.setParam('ccdread',"rdnoise")
    iraf.datapars.setParam('exposure',"exptime")
    iraf.datapars.setParam('airmass',"airmass")
    iraf.datapars.setParam('filter',"filter")
    iraf.datapars.setParam('obstime',"time-obs")
    # iraf.datapars.setParam('obstime',"date-obs")
    iraf.datapars.setParam('sigma',"INDEF")
    iraf.photpars.setParam('zmag',0.)
    iraf.centerpars.setParam('cbox',9.)
    iraf.centerpars.setParam('maxshift',3.)
    iraf.fitskypars.setParam('salgorithm',"median")
    iraf.fitskypars.setParam('dannulus',10.)
    
    # clean up the indefs so we can actually do stats, but reassign them to 99999 so we don't lose track of things
    # keep a separate list without them to do the median (we need floats)
    indefs = np.where(fwhm=='INDEF')
    good = np.where(fwhm!='INDEF')
    fwhm[indefs] = 99.999
    fwhm = fwhm.astype(float)
    fwhm_good = fwhm[good].astype(float)

    indefs = np.where(peak=='INDEF')
    peak[indefs] = -999.999
    peak = peak.astype(float)
    peak_good = peak[good].astype(float)
    
    if not os.path.isfile(image[0:-5]+'.txdump'):
        # findavgfwhm = sewpy.SEW(
        #     params = ["X_IMAGE", "Y_IMAGE", "FWHM_IMAGE", "FLAGS"],
        #     config = {"DETECT_THRESH":200.0},
        #     sexpath = "sex"
        # )
        # 
        # out = findavgfwhm(image)["table"]
        # 
        # fwhms = out['FWHM_IMAGE'] # This is an astropy table.
        # flags = out['FLAGS']
        
        # get a really rough estimate of the stellar FWHM in the image to set apertures
        
        # use the input fwhm measurement
        # ap1x = fwhm_est

        
        # xpos = datatable['X_IMAGE']
        # ypos = datatable['Y_IMAGE']
        # fwhm = datatable['FWHM_IMAGE']
        # flags = datatable['FLAGS']
        # idno = datatable['NUMBER']
        ap1x = np.median(fwhm_good) # only use isolated detections of stars, this is the 1x aperture
        # print ap1x
        ap2x = 2.0*ap1x
        
        # these = [ i for i,id in enumerate(idno) if (flags[i] == 0)]
        
        # with open(image[0:-5]+'.escut.pos','w+') as f:
        #     for j in range(len(xpos)):
        #         print >> f, xpos[j], ypos[j], fwhm[j], idno[j]
        
        iraf.datapars.setParam('fwhmpsf',ap1x)
        iraf.photpars.setParam('apertures',repr(ap1x)+', '+repr(ap2x))
        iraf.fitskypars.setParam('annulus',4.*ap1x)
        iraf.apphot.phot(image=image, coords=pos_file, output=image[0:-5]+'.phot')
        with open(image[0:-5]+'.txdump','w+') as txdump_out :
            iraf.ptools.txdump(textfiles=image[0:-5]+'.phot', fields="id,mag,merr,msky,stdev,rapert,xcen,ycen,ifilter,xairmass,image", expr='MAG[1] != INDEF && MERR[1] != INDEF && MAG[2] != INDEF && MERR[2] != INDEF', headers='no', Stdout=txdump_out)
            
    mag1x, mag2x = np.loadtxt(image[0:-5]+'.txdump', usecols=(1,2), unpack=True)
    iraf_id = np.loadtxt(image[0:-5]+'.txdump', usecols=(0,), dtype=int, unpack=True)
    # idno = np.loadtxt(image[0:-5]+'.escut.pos', usecols=(3,), dtype=int, unpack=True)
    xpos, ypos = np.loadtxt(pos_file, usecols=(0,1), unpack=True)
    
    keepIndex = iraf_id - 1
    
    xpos, ypos, fwhm, peak = xpos[keepIndex], ypos[keepIndex], fwhm[keepIndex], peak[keepIndex]
    
    # print idno.size, iraf_id.size, xpos.size
    
    diff = mag2x - mag1x
    
    diffCut = diff
    magCut = mag2x
    xCut = xpos#[good]
    yCut = ypos#[good]
    idCut = iraf_id
    fwhmCut = fwhm#_good
    peakCut = peak
    
    print peakCut.size, magCut.size, diffCut.size
    
    print diffCut.size, 0, np.median(diffCut), diffCut.std()
    nRemoved = 1   
    
    # plt.clf()
    # plt.scatter(peakCut, magCut, edgecolor='none')
    # plt.savefig('peaktest.pdf')
    
    plt.clf()
    # plt.hlines(bin_edges, -2, 1, colors='red', linestyle='dashed')
    plt.scatter(diff, mag2x, edgecolor='none', facecolor='black', s=4)
    # plt.scatter(diffCut, magCut, edgecolor='none', facecolor='blue', s=4)
    magdiff = zip(magCut.tolist(), diffCut.tolist(), peakCut.tolist(), idCut.tolist())
    dtype = [('mag',float), ('diff', float), ('peak', float), ('id', int)]
    magdiff = np.array(magdiff, dtype=dtype)
    
    magSort = np.sort(magdiff, order='peak')

    peakRange = (magSort['peak'] > 20000.0) & (magSort['peak'] < 40000.0)
    peakVal = np.median((magSort['diff'])[np.where(peakRange)])
    # peakVal = np.median(diffCut)
    print peakVal
    
    plt.scatter((magSort['diff'])[np.where(peakRange)], (magSort['mag'])[np.where(peakRange)], edgecolor='none', facecolor='blue', s=4)
    
    while nRemoved != 0:
        nBefore = diffCut.size
        diffCheck = np.where(abs(peakVal-diffCut) < 2.5*diffCut.std())#[i for i,d in enumerate(diff) if (-0.5 < d < 0.0)]
        # 
        diffCut = diffCut[diffCheck]
        nRemoved = nBefore - diffCut.size
        magCut = magCut[diffCheck]
        xCut = xCut[diffCheck]
        yCut = yCut[diffCheck]
        idCut = idCut[diffCheck]
        fwhmCut = fwhmCut[diffCheck]
        print diffCut.size, nRemoved, np.median(diffCut), diffCut.std()
        if 0.05 < diffCut.std() <0.06:
           nRemoved=0 
        # plt.fill_betweenx(bin_centers, bin_meds+3.0*bin_stds, bin_meds-3.0*bin_stds, facecolor='red', edgecolor='none', alpha=0.4, label='2x RMS sigma clipping region')
        
    # with open('escutSTD_i.pos','w+') as f:
    #     for i,blah in enumerate(xCut):
    #         print >> f, xCut[i], yCut[i], diffCut[i]
            
    bin_meds, bin_edges, binnumber = stats.binned_statistic(magCut, diffCut, statistic='median', bins=24, range=(magCut.min(),magCut.max()))
    bin_stds, bin_edges, binnumber = stats.binned_statistic(magCut, diffCut, statistic=np.std, bins=24, range=(magCut.min(),magCut.max()))
    bin_width = (bin_edges[1] - bin_edges[0])
    bin_centers = bin_edges[1:] - bin_width/2
    # print bin_meds, bin_stds
    
    left_edge = np.array(zip(peakVal-2.5*bin_stds, bin_centers))
    right_edge = np.flipud(np.array(zip(peakVal+2.5*bin_stds, bin_centers)))
    # print left_edge, right_edge
    verts = np.vstack((left_edge, right_edge))
    # print verts
    # verts = np.delete(verts, np.array([0,1,2,22,23,24,25,45,46,47]), axis=0)
    
    # print verts
    esRegion = Path(verts)
    sources = esRegion.contains_points(zip(diff,mag2x))
    # print sources
    
    with open('escutREG_i.pos','w+') as f:
        for i,blah in enumerate(xpos[sources]):
            print >> f, (xpos[sources])[i], (ypos[sources])[i], (diff[sources])[i]
    
    magCut = mag2x[sources]
    fwhmCut = fwhm[sources]   
    xCut = xpos[sources]
    yCut = ypos[sources] 
    diffCut = diff[sources]    
    
    # find the sources that are in the std method but not the region method
    # print idCut, idno[sources]
    # extrasSTD = np.setdiff1d(idno[sources], idCut)
    # print extrasSTD.size
    # print extrasSTD
    # with open('escutUNIQUE.pos','w+') as f:
    #     for i,blah in enumerate(extrasSTD):
    #         print >> f, xpos[blah-1], ypos[blah-1]
            
    # fwhmcheck = np.loadtxt('testfwhmREG.log', usecols=(10,), unpack=True)
    fwhmchk2 = np.where((magCut<-4) & (fwhmCut<90.0))
    print np.median(fwhmCut[fwhmchk2]), np.std(fwhmCut[fwhmchk2])
    fwchk = np.where(np.abs(fwhmCut-np.median(fwhmCut[fwhmchk2])) > 10.0*np.std(fwhmCut[fwhmchk2]))
    drop = np.abs(fwhmCut-np.median(fwhmCut[fwhmchk2])) > 10.0*np.std(fwhmCut[fwhmchk2])
    # fwmag = mag2x[sources]
    
    with open('escutVBAD_i.pos','w+') as f:
        for i,blah in enumerate(xCut[fwchk]):
            print >> f, (xCut[fwchk])[i], (yCut[fwchk])[i]
            
    with open('escut_r.pos','w+') as f:
        for i,blah in enumerate(xCut):
            if not drop[i]:
                print >> f, xCut[i], yCut[i]
    
    with open('escut_g.pos','w+') as f:
        for i,blah in enumerate(xCut):
            if not drop[i]:
                print >> f, xCut[i], yCut[i]
    
    plt.fill_betweenx(bin_centers, peakVal+2.5*bin_stds, peakVal-2.5*bin_stds, facecolor='red', edgecolor='none', alpha=0.4, label='2x RMS sigma clipping region')
    plt.scatter(diffCut[fwchk], magCut[fwchk], edgecolor='none', facecolor='red', s=4)
    plt.ylim(0,-12)
    plt.xlabel('$m_{2x} - m_{1x}$')
    plt.ylabel('$m_{2x}$')
    plt.xlim(-2,1)
    plt.savefig('testmagiraf.pdf')
    
    plt.clf()
    plt.scatter(magCut, fwhmCut, edgecolor='none', facecolor='black')
    plt.scatter(magCut[fwchk], fwhmCut[fwchk], edgecolor='none', facecolor='red')
    plt.hlines([np.median(fwhmCut)], -12, 0, colors='red', linestyle='dashed')
    plt.hlines([np.median(fwhmCut)+fwhmCut.std(), np.median(fwhmCut)-fwhmCut.std()], -12, 0, colors='red', linestyle='dotted')
    plt.ylim(0,20)
    plt.xlim(-12,0)
    plt.ylabel('fwhm')
    plt.xlabel('$m_{2x}$')
    plt.savefig('fwhmcheck.pdf')
    
    
    return True

