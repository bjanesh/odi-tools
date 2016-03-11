import sys, os, glob, string
import numpy as np
import astropy as ast
import matplotlib.pyplot as plt
from pyraf import iraf
import odi_config as odi

def get_sdss_coords(img, ota, output='test.sdss'):
    formats = ['csv','xml','html']

    astro_url='http://skyserver.sdss3.org/public/en/tools/search/x_sql.aspx'
    public_url='http://skyserver.sdss3.org/public/en/tools/search/x_sql.aspx'

    default_url=public_url
    default_fmt='csv'

    hdulist = odi.fits.open(img)
    hdu = hdulist[ota]
    pvlist = hdu.header['PV*']
    for pv in pvlist:
        tpv = 'T'+pv
        hdu.header.rename_keyword(pv, tpv, force=False)
    xdim = hdu.header['NAXIS1']
    ydim = hdu.header['NAXIS2']

    if not os.path.isfile(output):
        # and find the image center
        xc = xdim/2.0
        yc = ydim/2.0

        # get the CD matrix keywords
        cd11 = hdu.header['CD1_1']
        cd22 = hdu.header['CD2_2']
        # try to load cd12 and cd21, if they don't exist, set them to zero
        try :
            cd12 = hdu.header['CD1_2']
        except:
            cd12 = 0.0
        try :
            cd21 = hdu.header['CD2_1']
        except:
            cd21 = 0.0

        # print xdim, ydim, cd12, cd21

        # Parse the WCS keywords in the primary HDU
        w = odi.WCS(hdu.header)

        # Some pixel coordinates of interest.
        pixcrd = np.array([[xc,yc]], np.float_)

        # Convert pixel coordinates to world coordinates
        # The second argument is "origin" -- in this case we're declaring we
        # have 1-based (Fortran-like) coordinates.
        world = w.wcs_pix2world(pixcrd, 1)
        # print(world)    
        rac = world[0][0]
        decc = world[0][1]
        # print xc, yc, rac, decc

        # get the biggest radius of the image in arcminutes
        pixscal1 = 3600*abs(cd11)
        pixscal2 = 3600*abs(cd22)
        xas = pixscal1 * xdim
        yas = pixscal2 * ydim
        xam = xas/60
        yam = yas/60
        #print(xam,yam)
        #radius for query: sqrt2 = 1.414
        sizeam = 1.414*(xam+yam)/4
        print sizeam

        qry = "select O.ra, O.dec, O.psfMag_u, O.psfMagErr_u, O.psfMag_g, \nO.psfMagErr_g, O.psfMag_r, O.psfMagErr_r, O.psfMag_i, \nO.psfMagErr_i, O.psfMag_z, O.psfMagErr_z, O.probPSF \nfrom \ndbo.fGetNearbyObjEq("+repr(rac)+","+repr(decc)+","+repr(sizeam)+") \nas N inner join PhotoObjAll as O on O.objID = N.objID order by N.distance"

        #print 'with query\n-->', qry
        print 'fetching SDSS sources around',rac,decc,'with radius',sizeam,'arcmin'
        url = default_url
        fmt = default_fmt
        writefirst = 1
        verbose = 0

        ofp = open(output,'w+')
        if verbose:
            odi.write_header(ofp,'#',url,qry)
        file_ = odi.query(qry,url,fmt)
        # Output line by line (in case it's big)
        line = file_.readline()
        if line.startswith("ERROR"): # SQL Statement Error -> stderr
            ofp = sys.stderr
        if writefirst:
            ofp.write(string.rstrip(line)+os.linesep)
        line = file_.readline()
        while line:
            ofp.write(string.rstrip(line)+os.linesep)
            line = file_.readline()
        ofp.close()
    else:
        print 'SDSS sources already fetched!'
    return xdim, ydim
    hdulist.close()
    
def refetch_sdss_coords(img, ota, gapmask, gmaglim=19.,offline = False,source='sdss'):
    
    image = odi.reprojpath+'reproj_'+ota+'.'+str(img[16:])
    outcoords = odi.coordspath+'reproj_'+ota+'.'+str(img[16:-5])+'.sdss'
    hdulist = odi.fits.open(image)
    pvlist = hdulist[0].header['PV*']
    for pv in pvlist:
	tpv = 'T'+pv
	hdulist[0].header.rename_keyword(pv, tpv, force=False)
    xdim = hdulist[0].header['NAXIS1']
    ydim = hdulist[0].header['NAXIS2']  
  
    if offline == False:
	formats = ['csv','xml','html']
	astro_url='http://skyserver.sdss3.org/public/en/tools/search/x_sql.aspx'
	public_url='http://skyserver.sdss3.org/public/en/tools/search/x_sql.aspx'
	default_url=public_url
	default_fmt='csv'
	if not os.path.isfile(outcoords):
	    # and find the image center
	    xc = xdim/2.0
	    yc = ydim/2.0
	
	    # get the CD matrix keywords
	    cd11 = hdulist[0].header['CD1_1']
	    cd22 = hdulist[0].header['CD2_2']
	    # try to load cd12 and cd21, if they don't exist, set them to zero
	    try :
		cd12 = hdulist[0].header['CD1_2']
	    except:
		cd12 = 0.0
	    try :
		cd21 = hdulist[0].header['CD2_1']
	    except:
		cd21 = 0.0
	    
	    w = odi.WCS(hdulist[0].header)
	
	    # Some pixel coordinates of interest.
	    pixcrd = np.array([[xc,yc]], np.float_)
	    
	    # Convert pixel coordinates to world coordinates
	    # The second argument is "origin" -- in this case we're declaring we
	    # have 1-based (Fortran-like) coordinates.
	    world = w.wcs_pix2world(pixcrd, 1)
	    # print(world)    
	    rac = world[0][0]
	    decc = world[0][1]
	    # print xc, yc, rac, decc

	    # get the biggest radius of the image in arcminutes
	    pixscal1 = 3600*abs(cd11)
	    pixscal2 = 3600*abs(cd22)
	    xas = pixscal1 * xdim
	    yas = pixscal2 * ydim
	    xam = xas/60
	    yam = yas/60
	    #print(xam,yam)
	    #radius for query: sqrt2 = 1.414
	    sizeam = 1.414*(xam+yam)/4
	    print sizeam

	    qry = "select O.ra, O.dec, O.psfMag_u, O.psfMagErr_u, O.psfMag_g, \nO.psfMagErr_g, O.psfMag_r, O.psfMagErr_r, O.psfMag_i, \nO.psfMagErr_i, O.psfMag_z, O.psfMagErr_z, O.probPSF \nfrom \ndbo.fGetNearbyObjEq("+repr(rac)+","+repr(decc)+","+repr(sizeam)+") \nas N inner join PhotoObjAll as O on O.objID = N.objID order by N.distance"

	    #print 'with query\n-->', qry
	    print 'fetching SDSS sources around',rac,decc,'with radius',sizeam,'arcmin'
	    url = default_url
	    fmt = default_fmt
	    writefirst = 1
	    verbose = 0

	    ofp = open(outcoords,'w+')
	    if verbose:
		odi.write_header(ofp,'#',url,qry)
	    file_ = odi.query(qry,url,fmt)
	    # Output line by line (in case it's big)
	    line = file_.readline()
	    if line.startswith("ERROR"): # SQL Statement Error -> stderr
		ofp = sys.stderr
	    if writefirst:
		ofp.write(string.rstrip(line)+os.linesep)
	    line = file_.readline()
	    while line:
		ofp.write(string.rstrip(line)+os.linesep)
		line = file_.readline()
	    ofp.close()
	else:
	    print 'SDSS sources already fetched!'

	ras,decs,psfMag_u,psfMagErr_u,psfMag_g,psfMagErr_g,psfMag_r,psfMagErr_r,psfMag_i,psfMagErr_i,psfMag_z,psfMagErr_z = np.loadtxt(outcoords,usecols=(0,1,2,3,4,5,6,7,8,9,10,11), unpack=True, delimiter=',', skiprows=2)
	probPSF = np.loadtxt(outcoords, usecols=(12,), dtype=int, unpack=True, delimiter=',', skiprows=2)
	w = odi.WCS(hdulist[0].header)
	with open(odi.coordspath+'reproj_'+ota+'.'+str(img[16:-5])+'.sdssxy', 'w+') as fxy:
	    j=0
	    # k=0
	    for i,c in enumerate(ras):
		coords2 = [[ras[i],decs[i]]]
		pixcrd2 = w.wcs_world2pix(coords2, 1)
		if psfMag_g[i]<gmaglim and probPSF[i]==1:
		    if  100.0 <= pixcrd2[0][0] < xdim-100.0 and 100.0 <= pixcrd2[0][1] < ydim-100.0:
			# make an image cutout of the gap mask
			x, y = int(round(pixcrd2[0][0])), int(round(pixcrd2[0][1]))
			cutout = gapmask[y-30:y+30,x-30:x+30]
			# print cutout.flatten()
			# k+=1
			# print k,cutout.astype(bool).any()
		    if not (cutout.astype(bool)).any():
			# j+=1
			# print pixcrd2[0][0], pixcrd2[0][1],x-30,x+30,y-30,y+30
			# plt.imshow(cutout)
			# plt.show()
			print >> fxy, pixcrd2[0][0], pixcrd2[0][1], ras[i],decs[i],psfMag_u[i],psfMagErr_u[i],psfMag_g[i],psfMagErr_g[i],psfMag_r[i],psfMagErr_r[i],psfMag_i[i],psfMagErr_i[i],psfMag_z[i],psfMagErr_z[i]


    if offline == True:
	if source == 'sdss':
	    outcoords =  odi.sdsspath+'offline_'+ota+'.'+str(img[16:-5])+'.sdss'
	    ras,decs,psfMag_u,psfMagErr_u,psfMag_g,psfMagErr_g,psfMag_r,psfMagErr_r,psfMag_i,psfMagErr_i,psfMag_z,psfMagErr_z = np.loadtxt(outcoords,usecols=(0,1,2,3,4,5,6,7,8,9,10,11), unpack=True, delimiter=',', skiprows=1)

	if source == 'twomass':
	    outcoords =  odi.twomasspath+'offline_'+ota+'.'+str(img[16:-5])+'.mass'
	    ras,decs = np.loadtxt(outcoords,usecols=(2,3), unpack=True, delimiter=',', skiprows=1)
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
	print 'Using Ra and Dec from:',outcoords, 'for reproject'
	w = odi.WCS(hdulist[0].header)
	with open(odi.coordspath+'reproj_'+ota+'.'+str(img[16:-5])+'.sdssxy', 'w+') as fxy:
	    for i,c in enumerate(ras):
		coords2 = [[ras[i],decs[i]]]
		pixcrd2 = w.wcs_world2pix(coords2, 1)
		if psfMag_g[i]<gmaglim:
		    if  100.0 <= pixcrd2[0][0] < xdim-100.0 and 100.0 <= pixcrd2[0][1] < ydim-100.0:
			# make an image cutout of the gap mask
			x, y = int(round(pixcrd2[0][0])), int(round(pixcrd2[0][1]))
			cutout = gapmask[y-30:y+30,x-30:x+30]
			if not (cutout.astype(bool)).any():
			    print >> fxy, pixcrd2[0][0], pixcrd2[0][1], ras[i],decs[i],psfMag_u[i],psfMagErr_u[i],psfMag_g[i],psfMagErr_g[i],psfMag_r[i],psfMagErr_r[i],psfMag_i[i],psfMagErr_i[i],psfMag_z[i],psfMagErr_z[i]


    hdulist.close()
    
def repoxy_offline(img, ota, gapmask, gmaglim=19.,source='sdss'):
    image = odi.reprojpath+'reproj_'+ota+'.'+str(img[16:])
    hdulist = odi.fits.open(image)
    pvlist = hdulist[0].header['PV*']
    for pv in pvlist:
	tpv = 'T'+pv
	hdulist[0].header.rename_keyword(pv, tpv, force=False)
    xdim = hdulist[0].header['NAXIS1']
    ydim = hdulist[0].header['NAXIS2']
    if source == 'sdss':
	outcoords =  odi.sdsspath+'offline_'+ota+'.'+str(img[16:-5])+'.sdss'
	ras,decs,psfMag_u,psfMagErr_u,psfMag_g,psfMagErr_g,psfMag_r,psfMagErr_r,psfMag_i,psfMagErr_i,psfMag_z,psfMagErr_z = np.loadtxt(outcoords,usecols=(0,1,2,3,4,5,6,7,8,9,10,11), unpack=True, delimiter=',', skiprows=1)
	outputxy =  odi.coordspath+'reproj_'+ota+'.'+str(img[16:-5])+'.sdssxy'
    if source == 'twomass':
	outcoords =  odi.twomasspath+'offline_'+ota+'.'+str(img[16:-5])+'.mass'
	outputxy =  odi.coordspath+'reproj_'+ota+'.'+str(img[16:-5])+'.massxy'
	ras,decs = np.loadtxt(outcoords,usecols=(2,3), unpack=True, delimiter=',', skiprows=1)
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
    print 'Using Ra and Dec from:',outcoords, 'for reproject'
    w = odi.WCS(hdulist[0].header)
    with open(outputxy, 'w+') as fxy:
	    for i,c in enumerate(ras):
		coords2 = [[ras[i],decs[i]]]
		pixcrd2 = w.wcs_world2pix(coords2, 1)
		if psfMag_g[i]<gmaglim:
		    if  100.0 <= pixcrd2[0][0] < xdim-100.0 and 100.0 <= pixcrd2[0][1] < ydim-100.0:
			# make an image cutout of the gap mask
			x, y = int(round(pixcrd2[0][0])), int(round(pixcrd2[0][1]))
			cutout = gapmask[y-30:y+30,x-30:x+30]
			if not (cutout.astype(bool)).any():
			    print >> fxy, pixcrd2[0][0], pixcrd2[0][1], ras[i],decs[i],psfMag_u[i],psfMagErr_u[i],psfMag_g[i],psfMagErr_g[i],psfMag_r[i],psfMagErr_r[i],psfMag_i[i],psfMagErr_i[i],psfMag_z[i],psfMagErr_z[i]


    hdulist.close()

def sdss_coords_full(img, gmaglim=19.):
    formats = ['csv','xml','html']

    astro_url='http://skyserver.sdss3.org/public/en/tools/search/x_sql.aspx'
    public_url='http://skyserver.sdss3.org/public/en/tools/search/x_sql.aspx'

    default_url=public_url
    default_fmt='csv'

    image = img
    outcoords = img[:-5]+'.sdss'

    hdulist = odi.fits.open(image)
    data = hdulist[0].data
    # hdu = hdulist[ota]
    pvlist = hdulist[0].header['PV*']
    for pv in pvlist:
      tpv = 'T'+pv
      hdulist[0].header.rename_keyword(pv, tpv, force=False)
    xdim = hdulist[0].header['NAXIS1']
    ydim = hdulist[0].header['NAXIS2']

    if not os.path.isfile(outcoords):
      # and find the image center
      xc = xdim/2.0
      yc = ydim/2.0

      # get the CD matrix keywords
      cd11 = hdulist[0].header['CD1_1']
      cd22 = hdulist[0].header['CD2_2']
      # try to load cd12 and cd21, if they don't exist, set them to zero
      try :
          cd12 = hdulist[0].header['CD1_2']
      except:
          cd12 = 0.0
      try :
          cd21 = hdulist[0].header['CD2_1']
      except:
          cd21 = 0.0

      # print xdim, ydim, cd12, cd21

      # Parse the WCS keywords in the primary HDU
      w = odi.WCS(hdulist[0].header)

      # Some pixel coordinates of interest.
      pixcrd = np.array([[xc,yc]], np.float_)

      # Convert pixel coordinates to world coordinates
      # The second argument is "origin" -- in this case we're declaring we
      # have 1-based (Fortran-like) coordinates.
      world = w.wcs_pix2world(pixcrd, 1)
      # print(world)    
      rac = world[0][0]
      decc = world[0][1]
      # print xc, yc, rac, decc

      # get the biggest radius of the image in arcminutes
      pixscal1 = 3600*abs(cd11)
      pixscal2 = 3600*abs(cd22)
      xas = pixscal1 * xdim
      yas = pixscal2 * ydim
      xam = xas/60
      yam = yas/60
      #print(xam,yam)
      #radius for query: sqrt2 = 1.414
      sizeam = 1.414*(xam+yam)/4
      print sizeam

      qry = "select O.ra, O.dec, O.psfMag_u, O.psfMagErr_u, O.psfMag_g, \nO.psfMagErr_g, O.psfMag_r, O.psfMagErr_r, O.psfMag_i, \nO.psfMagErr_i, O.psfMag_z, O.psfMagErr_z, O.probPSF \nfrom \ndbo.fGetNearbyObjEq("+repr(rac)+","+repr(decc)+","+repr(sizeam)+") \nas N inner join PhotoObjAll as O on O.objID = N.objID order by N.distance"

      #print 'with query\n-->', qry
      print 'fetching SDSS sources around',rac,decc,'with radius',sizeam,'arcmin'
      url = default_url
      fmt = default_fmt
      writefirst = 1
      verbose = 0

      ofp = open(outcoords,'w+')
      if verbose:
          odi.write_header(ofp,'#',url,qry)
      file_ = odi.query(qry,url,fmt)
      # Output line by line (in case it's big)
      line = file_.readline()
      if line.startswith("ERROR"): # SQL Statement Error -> stderr
          ofp = sys.stderr
      if writefirst:
          ofp.write(string.rstrip(line)+os.linesep)
      line = file_.readline()
      while line:
          ofp.write(string.rstrip(line)+os.linesep)
          line = file_.readline()
      ofp.close()
    else:
      print 'SDSS sources already fetched!'

    ras,decs,psfMag_u,psfMagErr_u,psfMag_g,psfMagErr_g,psfMag_r,psfMagErr_r,psfMag_i,psfMagErr_i,psfMag_z,psfMagErr_z = np.loadtxt(outcoords,usecols=(0,1,2,3,4,5,6,7,8,9,10,11), unpack=True, delimiter=',', skiprows=2)
    probPSF = np.loadtxt(outcoords, usecols=(12,), dtype=int, unpack=True, delimiter=',', skiprows=2)
    # print ras, decs

    w = odi.WCS(hdulist[0].header)
    with open(img[:-5]+'.wcs.coo','w+') as f:
        with open(img[:-5]+'.sdssxy', 'w+') as fxy:
            for i,c in enumerate(ras):
                coords2 = [[ras[i],decs[i]]]
                pixcrd2 = w.wcs_world2pix(coords2, 1)
                if psfMag_g[i]<gmaglim and probPSF[i]==1 and 100.0 <= pixcrd2[0][0] < xdim-100.0 and 100.0 <= pixcrd2[0][1] < ydim-100.0:
                    r, d = odi.deg_to_sex(ras[i], decs[i])
                    x, y = int(round(pixcrd2[0][0])), int(round(pixcrd2[0][1]))
                    cutout = data[y-30:y+30,x-30:x+30]
                    # print cutout.flatten()
                    # k+=1
                    # print i,(cutout<-900).any()
                    if not (cutout<-900).any():
                      print >> f, r, d, psfMag_g[i]
                      print >> fxy, pixcrd2[0][0], pixcrd2[0][1], ras[i],decs[i],psfMag_u[i],psfMagErr_u[i],psfMag_g[i],psfMagErr_g[i],psfMag_r[i],psfMagErr_r[i],psfMag_i[i],psfMagErr_i[i],psfMag_z[i],psfMagErr_z[i]
    hdulist.close()    

def get_sdss_coords_offline(img, ota, output='test.sdss'):
    hdulist = odi.fits.open(img)
    hdu = hdulist[ota]
    
    pvlist = hdu.header['PV*']
    for pv in pvlist:
        tpv = 'T'+pv
        hdu.header.rename_keyword(pv, tpv, force=False)
    xdim = hdu.header['NAXIS1']
    ydim = hdu.header['NAXIS2']
    
    sdss_cat_img = hdulist['CAT.PHOTCALIB']
    sdss_cat_img_df = pd.DataFrame.from_dict(sdss_cat_img.data)
    hdulist1.close()
    
    ota = float(ota.strip('OTA.SCI'))
    ota_matches_df = sdss_cat_img_df.iloc[np.where(sdss_cat_img_df['ODI_OTA'] == ota)]

    needed_columns = ['SDSS_RA','SDSS_DEC','SDSS_MAG_U',
		      'SDSS_ERR_U', u'SDSS_MAG_G', u'SDSS_ERR_G', u'SDSS_MAG_R',
		      'SDSS_ERR_R', u'SDSS_MAG_I', u'SDSS_ERR_I', u'SDSS_MAG_Z',
		      'SDSS_ERR_Z','ODI_OTA','ODI_OTA']

    output_df = ota_matches_df[needed_columns]
    output_df.to_csv(output,index=False)
    return xdim, ydim

def refetch_sdss_coords_offline(img, ota, gapmask, gmaglim=19.):
    image = odi.reprojpath+'reproj_'+ota+'.'+str(img[16:])
    outcoords = odi.coordspath+'reproj_'+ota+'.'+str(img[16:-5])+'.sdss'

    hdulist = odi.fits.open(image)
    pvlist = hdulist[0].header['PV*']
    for pv in pvlist:
      tpv = 'T'+pv
      hdulist[0].header.rename_keyword(pv, tpv, force=False)
    xdim = hdulist[0].header['NAXIS1']
    ydim = hdulist[0].header['NAXIS2']

    if not os.path.isfile(outcoords):
      xc = xdim/2.0
      yc = ydim/2.0

      # get the CD matrix keywords
      cd11 = hdulist[0].header['CD1_1']
      cd22 = hdulist[0].header['CD2_2']
      # try to load cd12 and cd21, if they don't exist, set them to zero
      try :
          cd12 = hdulist[0].header['CD1_2']
      except:
          cd12 = 0.0
      try :
          cd21 = hdulist[0].header['CD2_1']
      except:
          cd21 = 0.0

      # print xdim, ydim, cd12, cd21

      # Parse the WCS keywords in the primary HDU
      w = odi.WCS(hdulist[0].header)

      # Some pixel coordinates of interest.
      pixcrd = np.array([[xc,yc]], np.float_)

      # Convert pixel coordinates to world coordinates
      # The second argument is "origin" -- in this case we're declaring we
      # have 1-based (Fortran-like) coordinates.
      world = w.wcs_pix2world(pixcrd, 1)
      # print(world)    
      rac = world[0][0]
      decc = world[0][1]
      # print xc, yc, rac, decc

      # get the biggest radius of the image in arcminutes
      pixscal1 = 3600*abs(cd11)
      pixscal2 = 3600*abs(cd22)
      xas = pixscal1 * xdim
      yas = pixscal2 * ydim
      xam = xas/60
      yam = yas/60
      #print(xam,yam)
      #radius for query: sqrt2 = 1.414
      sizeam = 1.414*(xam+yam)/4
      print sizeam



    ras,decs,psfMag_u,psfMagErr_u,psfMag_g,psfMagErr_g,psfMag_r,psfMagErr_r,psfMag_i,psfMagErr_i,psfMag_z,psfMagErr_z = np.loadtxt(outcoords,usecols=(0,1,2,3,4,5,6,7,8,9,10,11), unpack=True, delimiter=',', skiprows=2)
    probPSF = np.loadtxt(outcoords, usecols=(12,), dtype=int, unpack=True, delimiter=',', skiprows=2)
    # print ras, decs

    w = odi.WCS(hdulist[0].header)

    with open(odi.coordspath+'reproj_'+ota+'.'+str(img[16:-5])+'.sdssxy', 'w+') as fxy:
        j=0
        # k=0
        for i,c in enumerate(ras):
          coords2 = [[ras[i],decs[i]]]
          pixcrd2 = w.wcs_world2pix(coords2, 1)
          if psfMag_g[i]<gmaglim and probPSF[i]==1:
            if  100.0 <= pixcrd2[0][0] < xdim-100.0 and 100.0 <= pixcrd2[0][1] < ydim-100.0:
                # make an image cutout of the gap mask
                x, y = int(round(pixcrd2[0][0])), int(round(pixcrd2[0][1]))
                cutout = gapmask[y-30:y+30,x-30:x+30]
                # print cutout.flatten()
                # k+=1
                # print k,cutout.astype(bool).any()
                if not (cutout.astype(bool)).any():
                    # j+=1
                    # print pixcrd2[0][0], pixcrd2[0][1],x-30,x+30,y-30,y+30
                    # plt.imshow(cutout)
                    # plt.show()
                    print >> fxy, pixcrd2[0][0], pixcrd2[0][1], ras[i],decs[i],psfMag_u[i],psfMagErr_u[i],psfMag_g[i],psfMagErr_g[i],psfMag_r[i],psfMagErr_r[i],psfMag_i[i],psfMagErr_i[i],psfMag_z[i],psfMagErr_z[i]
        # print j
    hdulist.close()
    

if __name__ == '__main__':
    get_sdss_coords()

