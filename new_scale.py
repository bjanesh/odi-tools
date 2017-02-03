import sys, os, glob, string
import pandas as pd
from collections import OrderedDict
import shutil
import numpy as np
from pyraf import iraf
import odi_config as odi


def img_sdss(ref_img,img,source='sdss'):
    if source == 'sdss':
	hdulist1 = odi.fits.open(ref_img)
	sdss_cat_ref_img = hdulist1['CAT.PHOTCALIB']
	sdss_cat_ref_img_df = pd.DataFrame.from_dict(sdss_cat_ref_img.data)
	hdulist1.close()
    
	hdulist2 = odi.fits.open(img)
	sdss_cat_img = hdulist2['CAT.PHOTCALIB']
	sdss_cat_img_df = pd.DataFrame.from_dict(sdss_cat_img.data)
	hdulist2.close()
    
	img_match = pd.merge(sdss_cat_ref_img_df,sdss_cat_img_df ,on = ['SDSS_RA','SDSS_DEC'],how='inner',suffixes=['_ref', '_img'])
    
	ota_matches = {}
	for key in odi.OTA_dictionary:
	    ota = float(odi.OTA_dictionary[key].strip('OTA.SCI'))
	    ota_matches['sdss_'+str(ota).strip('.0')] = img_match.iloc[np.where((img_match['ODI_OTA_img'] == ota) & (img_match['ODI_OTA_ref'] == ota))]

	needed_columns = ['SDSS_RA','SDSS_DEC','SDSS_MAG_U_img',
			'SDSS_ERR_U_img', u'SDSS_MAG_G_img', u'SDSS_ERR_G_img', u'SDSS_MAG_R_img',
			'SDSS_ERR_R_img', u'SDSS_MAG_I_img', u'SDSS_ERR_I_img', u'SDSS_MAG_Z_img',
			'SDSS_ERR_Z_img','ODI_OTA_ref','ODI_OTA_img']
	output = img_match[needed_columns]
	regtest = img_match[['SDSS_RA','SDSS_DEC']]
	regtest.to_csv(odi.matchpath+img.split('.')[1][0] + '-'+ref_img.split('.')[1][0]+'.reg',index=False)
	match_file = odi.matchpath+img.split('.')[1][0] + '-'+ref_img.split('.')[1][0]+'.match'
	output.to_csv(match_file,index=False)
	return img_match,ota_matches
    if source == 'twomass':
	hdulist1 = odi.fits.open(ref_img)
	mass_cat_ref = hdulist1['CAT.ODI+2MASS']
	ref_ra =  mass_cat_ref.data['TWOMASS_RA']
	ref_dec = mass_cat_ref.data['TWOMASS_DEC']
	ref_ota = mass_cat_ref.data['OTA']
	
	ref_values = {'TWOMASS_RA':ref_ra,'TWOMASS_DEC':ref_dec,'OTA':ref_ota}
    
	mass_cat_ref_df = pd.DataFrame.from_dict(ref_values)
	new_ref = mass_cat_ref_df.apply(lambda x: x.values.byteswap().newbyteorder())
	hdulist1.close()
    
	hdulist2 = odi.fits.open(img)
	mass_cat_img = hdulist2['CAT.ODI+2MASS']
	img_ra =  mass_cat_img.data['TWOMASS_RA']
	img_dec = mass_cat_img.data['TWOMASS_DEC']
	img_ota = mass_cat_img.data['OTA']
	img_values = {'TWOMASS_RA':img_ra,'TWOMASS_DEC':img_dec,'OTA':img_ota}
    
	mass_cat_img_df = pd.DataFrame.from_dict(img_values)
	new_img = mass_cat_img_df.apply(lambda x: x.values.byteswap().newbyteorder())
	hdulist2.close()
    
    
	img_match = pd.merge(new_ref,new_img ,on = ['TWOMASS_RA','TWOMASS_DEC'],how='inner',suffixes=['_ref', '_img'])
	needed_columns = ['TWOMASS_RA','TWOMASS_DEC','OTA_ref','OTA_img']
	output = img_match[needed_columns]

	match_file = odi.matchpath+img.split('.')[1][0] + '-'+ref_img.split('.')[1][0]+'.mass.match'
	output.to_csv(match_file,index=False)
	return img_match

def img_twomass(ref_img,img):
    hdulist1 = odi.fits.open(ref_img)
    mass_cat_ref = hdulist1['CAT.ODI+2MASS']
    ref_ra =  mass_cat_ref.data['TWOMASS_RA']
    ref_dec = mass_cat_ref.data['TWOMASS_DEC']
    ref_ota = mass_cat_ref.data['OTA']
    
    ref_values = {'TWOMASS_RA':ref_ra,'TWOMASS_DEC':ref_dec,'OTA':ref_ota}
    
    mass_cat_ref_df = pd.DataFrame.from_dict(ref_values)
    new_ref = mass_cat_ref_df.apply(lambda x: x.values.byteswap().newbyteorder())
    #print mass_cat_ref_df

    hdulist1.close()
    
    hdulist2 = odi.fits.open(img)
    mass_cat_img = hdulist2['CAT.ODI+2MASS']
    img_ra =  mass_cat_img.data['TWOMASS_RA']
    img_dec = mass_cat_img.data['TWOMASS_DEC']
    img_ota = mass_cat_img.data['OTA']
    
    psfMag_u       = np.ones(len(img_ra))
    psfMagErr_u    = np.ones(len(img_ra))
    psfMag_g       = np.ones(len(img_ra))
    psfMagErr_g    = np.ones(len(img_ra))
    psfMag_r       = np.ones(len(img_ra))
    psfMagErr_r    = np.ones(len(img_ra))
    psfMag_i       = np.ones(len(img_ra))
    psfMagErr_i    = np.ones(len(img_ra))
    psfMag_z       = np.ones(len(img_ra))
    psfMagErr_z    = np.ones(len(img_ra))
    
    img_values = {'TWOMASS_RA':img_ra,'TWOMASS_DEC':img_dec,'OTA':img_ota}
    
    mass_cat_img_df = pd.DataFrame.from_dict(img_values)
    new_img = mass_cat_img_df.apply(lambda x: x.values.byteswap().newbyteorder())
    #comp_img_df.dtypes
    hdulist2.close()
    
    
    img_match = pd.merge(new_ref,new_img ,on = ['TWOMASS_RA','TWOMASS_DEC'],how='inner',suffixes=['_ref', '_img'])
    needed_columns = ['TWOMASS_RA','TWOMASS_DEC','OTA_ref','OTA_img']
    output = img_match[needed_columns]

    match_file = odi.matchpath+img.split('.')[1][0] + '-'+ref_img.split('.')[1][0]+'.mass.match'
    output.to_csv(match_file,index=False)
    return img_match


def list_wcs_coords_match(img,ref_img, ota, gapmask, gmaglim=23.):
    image = odi.reprojpath+'reproj_'+ota+'.'+img.stem()
    outcoords = odi.matchpath+img.split('.')[1][0] + '-'+ref_img.split('.')[1][0]+'.match'
    
    hdulist = ast.io.odi.fits.open(image)
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
    ras,decs,psfMag_u,psfMagErr_u,psfMag_g,psfMagErr_g,psfMag_r,psfMagErr_r,psfMag_i,psfMagErr_i,psfMag_z,psfMagErr_z = np.loadtxt(outcoords,usecols=(0,1,2,3,4,5,6,7,8,9,10,11), unpack=True, delimiter=',', skiprows=1)
    
    w = odi.WCS(hdulist[0].header)
    
    with open(odi.matchpath+'reproj_'+ota+'.'+img.base()+'.sdssxy', 'w+') as fxy:
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

def match_xy(img,ref_img,ota,gapmask, gmaglim=23., source ='sdss'):
    suffix = '.'+source.strip('two')+'.match'
    outcoords = odi.matchpath+img.split('.')[1][0] + '-'+ref_img.split('.')[1][0]+suffix
    if source == 'sdss':
	ras,decs,psfMag_u,psfMagErr_u,psfMag_g,psfMagErr_g,psfMag_r,psfMagErr_r,psfMag_i,psfMagErr_i,psfMag_z,psfMagErr_z,ODI_OTA_ref,ODI_OTA_img = np.loadtxt(outcoords,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13), unpack=True, delimiter=',', skiprows=1)
	outputxy  = odi.matchpath+'reproj_'+ota+'.'+img.base()+'.sdssxy'
    if source == 'twomass':
	outputxy  = odi.matchpath+'reproj_'+ota+'.'+img.base()+'.massxy'
	ras,decs,ODI_OTA_ref,ODI_OTA_img = np.loadtxt(outcoords,usecols=(0,1,2,3), unpack=True, delimiter=',', skiprows=1)
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
    stars_on_ota = np.where(ODI_OTA_img.astype(int) == int(ota.strip('OTA.SCI')))
    ras           =  ras[stars_on_ota[0]] 
    decs          =  decs[stars_on_ota[0]] 
    psfMag_u      =  psfMag_u[stars_on_ota[0]] 
    psfMagErr_u   =  psfMagErr_u[stars_on_ota[0]] 
    psfMag_g      =  psfMag_g[stars_on_ota[0]] 
    psfMagErr_g   =  psfMagErr_g[stars_on_ota[0]] 
    psfMag_r      =  psfMag_r[stars_on_ota[0]] 
    psfMagErr_r   =  psfMagErr_r[stars_on_ota[0]] 
    psfMag_i      =  psfMag_i[stars_on_ota[0]] 
    psfMagErr_i   =  psfMagErr_i[stars_on_ota[0]] 
    psfMag_z      =  psfMag_z[stars_on_ota[0]] 
    psfMagErr_z   =  psfMagErr_z[stars_on_ota[0]] 
    ODI_OTA_ref   =  ODI_OTA_ref[stars_on_ota[0]] 
    ODI_OTA_img   =  ODI_OTA_img[stars_on_ota[0]]

    image = odi.reprojpath+'reproj_'+ota+'.'+img.stem()
    hdulist = odi.fits.open(image)
    
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
    w = odi.WCS(hdulist[0].header)
    with open(outputxy, 'w+') as fxy:
	j=0
	# k=0
	for i,c in enumerate(ras):
	    coords2 = [[ras[i],decs[i]]]
	    pixcrd2 = w.wcs_world2pix(coords2, 1)
	    if psfMag_g[i]<gmaglim:
		if  100.0 <= pixcrd2[0][0] < xdim-100.0 and 100.0 <= pixcrd2[0][1] < ydim-100.0:
		    # make an image cutout of the gap mask
		    x, y = int(round(pixcrd2[0][0])), int(round(pixcrd2[0][1]))
		    cutout = gapmask[y-30:y+30,x-30:x+30]
		    if not (cutout.astype(bool)).any():
			print >> fxy, pixcrd2[0][0], pixcrd2[0][1], ras[i],decs[i],psfMag_u[i],psfMagErr_u[i],psfMag_g[i],psfMagErr_g[i],psfMag_r[i],psfMagErr_r[i],psfMag_i[i],psfMagErr_i[i],psfMag_z[i],psfMagErr_z[i],ODI_OTA_ref[i],ODI_OTA_img[i]
      # print j
    hdulist.close()
    if os.stat(outputxy).st_size != 0:
	match_status_1 = 'OK'
    else:
	match_status_1 = 'BAD'
    return match_status_1


def ref_xy(img,ref_img,ota,gmaglim=23.,source='sdss'):
    if source == 'sdss':
	suffix = '.sdssxy'
    if source == 'twomass':
	suffix = '.massxy'
    print 'reading',odi.matchpath+'reproj_'+ota+'.'+img.base()+suffix
    inputxy = odi.matchpath+'reproj_'+ota+'.'+img.base()+suffix
    pix1,pix2,ras,decs,psfMag_u,psfMagErr_u,psfMag_g,psfMagErr_g,psfMag_r,psfMagErr_r,psfMag_i,psfMagErr_i,psfMag_z,psfMagErr_z,ODI_OTA_ref,ODI_OTA_img = np.loadtxt(inputxy,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), unpack=True)
    
    if type(pix1) == np.float64:
	pix1,pix2,ras,decs,psfMag_u,psfMagErr_u,psfMag_g,psfMagErr_g,psfMag_r,psfMagErr_r,psfMag_i,psfMagErr_i,psfMag_z,psfMagErr_z,ODI_OTA_ref,ODI_OTA_img = np.loadtxt(inputxy,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), unpack=True).reshape((-1,1))
    
    for o in range(len(ODI_OTA_ref)):
	ota_ref = 'OTA'+str(int(ODI_OTA_ref[o]))+'.SCI'
	image = odi.reprojpath+'reproj_'+ota_ref+'.'+str(ref_ref_img.stem())
	hdulist = odi.fits.open(image)
	pvlist = hdulist[0].header['PV*']
	for pv in pvlist:
	    tpv = 'T'+pv
	    hdulist[0].header.rename_keyword(pv, tpv, force=False)
	xdim = hdulist[0].header['NAXIS1']
	ydim = hdulist[0].header['NAXIS2']
	w = odi.WCS(hdulist[0].header)
	
	
	gapmask = odi.get_gaps_rep(ref_img,ota_ref)
	
	
	ref_output = odi.matchpath+'reproj_'+ota+'.'+img.base()+suffix+'.ref'
	reg_output = odi.matchpath+'reproj_'+ota+'.'+img.base()+suffix+'.ref'+'.reg'
	with open(ref_output, 'w+') as fxy:
	    with open(reg_output,'w+') as rx:
		j=0
		for i,c in enumerate(ras):
		    coords2 = [[ras[i],decs[i]]]
		    pixcrd2 = w.wcs_world2pix(coords2, 1)
		    if psfMag_g[i]<gmaglim:
			if  100.0 <= pixcrd2[0][0] < xdim-100.0 and 100.0 <= pixcrd2[0][1] < ydim-100.0:
			    # make an image cutout of the gap mask
			    x, y = int(round(pixcrd2[0][0])), int(round(pixcrd2[0][1]))
			    cutout = gapmask[y-30:y+30,x-30:x+30]
			    if not (cutout.astype(bool)).any():
				print >> fxy, pixcrd2[0][0], pixcrd2[0][1], ras[i],decs[i],psfMag_u[i],psfMagErr_u[i],psfMag_g[i],psfMagErr_g[i],psfMag_r[i],psfMagErr_r[i],psfMag_i[i],psfMagErr_i[i],psfMag_z[i],psfMagErr_z[i],ODI_OTA_ref[i],ODI_OTA_img[i]
				print >> rx, pixcrd2[0][0], pixcrd2[0][1]
		    #if psfMag_g[i]<gmaglim:
			#if  100.0 <= pixcrd2[0][0] < xdim-100.0 and 100.0 <= pixcrd2[0][1] < ydim-100.0:
			    ## make an image cutout of the gap mask
			    #x, y = int(round(pixcrd2[0][0])), int(round(pixcrd2[0][1]))
			    #cutout = gapmask[y-30:y+30,x-30:x+30]
			    #if not (cutout.astype(bool)).any():
			#print >> fxy, pixcrd2[0][0], pixcrd2[0][1], ras[i],decs[i],psfMag_u[i],psfMagErr_u[i],psfMag_g[i],psfMagErr_g[i],psfMag_r[i],psfMagErr_r[i],psfMag_i[i],psfMagErr_i[i],psfMag_z[i],psfMagErr_z[i],ODI_OTA_ref[i],ODI_OTA_img[i]
	hdulist.close()
	if (os.stat(odi.matchpath+'reproj_'+ota+'.'+img.base()+suffix+'.ref').st_size != 0):
	    match_status = 'OK'
	else:
	    match_status = 'BAD'
    return match_status

def match_phot_coords(img,ref_img, ota,source='sdss'):
    if source == 'sdss':
	suffix = '.sdssxy'
    if source == 'twomass':
	suffix = '.massxy'
    ota_cat = pix1,pix2,ras,decs,psfMag_u,psfMagErr_u,psfMag_g,psfMagErr_g,psfMag_r,psfMagErr_r,psfMag_i,psfMagErr_i,psfMag_z,psfMagErr_z,ODI_OTA_ref,ODI_OTA_img  = np.loadtxt(odi.matchpath+'reproj_'+ota+'.'+img.base()+suffix,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), unpack=True)
    ref_cat = pix1_ref,pix2_ref,ras_ref,decs_ref,psfMag_u_ref,psfMagErr_u_ref,psfMag_g_ref,psfMagErr_g_ref,psfMag_r_ref,psfMagErr_r_ref,psfMag_i_ref,psfMagErr_i_ref,psfMag_z_ref,psfMagErr_z_ref,ODI_OTA_ref_ref,ODI_OTA_img_ref = np.loadtxt(odi.matchpath+'reproj_'+ota+'.'+img.base()+suffix+'.ref',usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), unpack=True)
    
    if type(ras) == np.float64:
	pix1,pix2,ras,decs,psfMag_u,psfMagErr_u,psfMag_g,psfMagErr_g,psfMag_r,psfMagErr_r,psfMag_i,psfMagErr_i,psfMag_z,psfMagErr_z,ODI_OTA_ref,ODI_OTA_img  = np.loadtxt(odi.matchpath+'reproj_'+ota+'.'+img.base()+suffix,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), unpack=True).reshape((-1,1))
    if type(ras_ref) == np.float64:
	pix1_ref,pix2_ref,ras_ref,decs_ref,psfMag_u_ref,psfMagErr_u_ref,psfMag_g_ref,psfMagErr_g_ref,psfMag_r_ref,psfMagErr_r_ref,psfMag_i_ref,psfMagErr_i_ref,psfMag_z_ref,psfMagErr_z_ref,ODI_OTA_ref_ref,ODI_OTA_img_ref = np.loadtxt(odi.matchpath+'reproj_'+ota+'.'+img.base()+suffix+'.ref',usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), unpack=True).reshape((-1,1))
    
    ota_cat_outfile = odi.matchpath+'reproj_'+ota+'.'+img.base()+suffix+'.coo'
    ref_cat_outfile = odi.matchpath+'reproj_'+ota+'.'+img.base()+suffix+'.ref.coo'
    with open(ota_cat_outfile, 'w+') as oxy:
	with open(ref_cat_outfile, 'w+') as rox:
	    for i in range(len(ras)):
		for j in range(len(ras_ref)):
		    if ((ras[i] == ras_ref[j]) and (decs[i] == decs_ref[j])):
			print >> oxy,pix1[i], pix2[i], ras[i],decs[i],psfMag_u[i],psfMagErr_u[i],psfMag_g[i],psfMagErr_g[i],psfMag_r[i] ,psfMagErr_r[i],psfMag_i[i],psfMagErr_i[i],psfMag_z[i],psfMagErr_z[i],ODI_OTA_ref[i],ODI_OTA_img[i]
			print >> rox,pix1_ref[j], pix2_ref[j], ras_ref[j],decs_ref[j],psfMag_u_ref[j],psfMagErr_u_ref[j],psfMag_g_ref[j],psfMagErr_g_ref[j],psfMag_r_ref[j],psfMagErr_r_ref[j],psfMag_i_ref[j],psfMagErr_i_ref[j],psfMag_z_ref[j],psfMagErr_z_ref[j],ODI_OTA_ref_ref[j],ODI_OTA_img_ref[j]
    return



def img_ref_fwhm(img,ref_img, ota,source = 'sdss'):
    if source == 'sdss':
	suffix = '.sdssxy'
    if source == 'twomass':
	suffix = '.massxy'    
    ota_coo_file = odi.matchpath+'reproj_'+ota+'.'+img.base()+suffix+'.coo'
    ref_coo_file = odi.matchpath+'reproj_'+ota+'.'+img.base()+suffix+'.ref.coo'
    ota_cat = pix1,pix2,ras,decs,psfMag_u,psfMagErr_u,psfMag_g,psfMagErr_g,psfMag_r,psfMagErr_r,psfMag_i,psfMagErr_i,psfMag_z,psfMagErr_z,ODI_OTA_ref,ODI_OTA_img  = np.loadtxt(ota_coo_file,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), unpack=True)
    ref_cat = pix1_ref,pix2_ref,ras_ref,decs_ref,psfMag_u_ref,psfMagErr_u_ref,psfMag_g_ref,psfMagErr_g_ref,psfMag_r_ref,psfMagErr_r_ref,psfMag_i_ref,psfMagErr_i_ref,psfMag_z_ref,psfMagErr_z_ref,ODI_OTA_ref_ref,ODI_OTA_img_ref = np.loadtxt(ref_coo_file,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), unpack=True)
    
    if type(ras) == np.float64:
	pix1,pix2,ras,decs,psfMag_u,psfMagErr_u,psfMag_g,psfMagErr_g,psfMag_r,psfMagErr_r,psfMag_i,psfMagErr_i,psfMag_z,psfMagErr_z,ODI_OTA_ref,ODI_OTA_img  = np.loadtxt(ota_coo_file,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), unpack=True).reshape((-1,1))
    if type(ras_ref) == np.float64:
	pix1_ref,pix2_ref,ras_ref,decs_ref,psfMag_u_ref,psfMagErr_u_ref,psfMag_g_ref,psfMagErr_g_ref,psfMag_r_ref,psfMagErr_r_ref,psfMag_i_ref,psfMagErr_i_ref,psfMag_z_ref,psfMagErr_z_ref,ODI_OTA_ref_ref,ODI_OTA_img_ref = np.loadtxt(ref_coo_file,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), unpack=True).reshape((-1,1))
    
    print len(set(ODI_OTA_ref_ref))
    
    radius=4.0 
    buff=7.0 
    width=5.0
    
    image = odi.bgsubpath+'bgsub_'+ota+'.'+img.stem()
    coords = odi.matchpath+'reproj_'+ota+'.'+img.base()+suffix+'.coo'
    outputfile = odi.matchpath+img.nofits()+'.'+ota+'.fwhm.log'

    print 'hello',image, coords
    
    iraf.tv.rimexam.setParam('radius',radius)
    iraf.tv.rimexam.setParam('buffer',buff)
    iraf.tv.rimexam.setParam('width',width)
    iraf.tv.rimexam.setParam('rplot',20.)
    iraf.tv.rimexam.setParam('center','yes')
    # fit a gaussian, rather than a moffat profile (it's more robust for faint sources)
    iraf.tv.rimexam.setParam('fittype','gaussian')
    iraf.tv.rimexam.setParam('iterati',1)
    
    if not os.path.isfile(outputfile):
        iraf.tv.imexamine(image, frame=10, logfile = outputfile, keeplog = 'yes', defkey = "a", nframes=0, imagecur = coords, wcs = "logical", use_display='no',  StdoutG='/dev/null',mode='h')
    outputfile_clean = open(outputfile.replace('.log','_clean.log'),"w")
    for line in open(outputfile,"r"):
      if not 'INDEF' in line:
        outputfile_clean.write(line)
      if 'INDEF' in line:
	outputfile_clean.write(line.replace('INDEF','999'))
    outputfile_clean.close()
    os.rename(outputfile.replace('.log','_clean.log'),outputfile)
    
    gfwhm_ota = np.loadtxt(outputfile, usecols=(10,), unpack=True)

    print 'median gwfhm in ota',ota+': ',np.median(gfwhm_ota),'pixels'# (determined via QR)'
    
    ota_ref = 'OTA'+str(int(ODI_OTA_ref_ref[0]))+'.SCI'
    image = odi.bgsubpath+'bgsub_'+ota_ref+'.'+str(ref_ref_img.stem())
    coords = odi.matchpath+'reproj_'+ota+'.'+img.base()+suffix+'.ref.coo'
    outputfile = odi.matchpath+img.nofits()+'.'+ota+'.fwhm.ref.log'
    
    print image, coords
    
    iraf.tv.rimexam.setParam('radius',radius)
    iraf.tv.rimexam.setParam('buffer',buff)
    iraf.tv.rimexam.setParam('width',width)
    iraf.tv.rimexam.setParam('rplot',20.)
    iraf.tv.rimexam.setParam('center','yes')
    # fit a gaussian, rather than a moffat profile (it's more robust for faint sources)
    iraf.tv.rimexam.setParam('fittype','gaussian')
    iraf.tv.rimexam.setParam('iterati',1)
    
    if not os.path.isfile(outputfile):
        iraf.tv.imexamine(image, frame=10, logfile = outputfile, keeplog = 'yes', defkey = "a", nframes=0, imagecur = coords, wcs = "logical", use_display='no',  StdoutG='/dev/null',mode='h')
    outputfile_clean = open(outputfile.replace('.log','_clean.log'),"w")
    for line in open(outputfile,"r"):
      if not 'INDEF' in line:
        outputfile_clean.write(line)
      if 'INDEF' in line:
	outputfile_clean.write(line.replace('INDEF','999'))
    outputfile_clean.close()
    os.rename(outputfile.replace('.log','_clean.log'),outputfile)
    
    gfwhm_ref = np.loadtxt(outputfile, usecols=(10,), unpack=True)

    print 'median gwfhm in ref ota',ota+': ',np.median(gfwhm_ref),'pixels'# (determined via QR)'
    
    return np.median(gfwhm_ota)



def new_getscale(img,ref_img,filter,verbose=False):
    sigThreshold = 0.006
    
    os.system("cat match/*.fwhm.log > match/img_fwhm.txt")
    
    ota_mag,ota_flux,ota_sky,ota_peak,ota_fwhm = np.loadtxt('match/img_fwhm.txt',usecols=(2,3,4,9,10),unpack=True)
    
    os.system("cat match/*.fwhm.ref.log > match/ref_fwhm.txt")
    
    ref_mag,ref_flux,ref_sky,ref_peak,ref_fwhm = np.loadtxt('match/ref_fwhm.txt',usecols=(2,3,4,9,10),unpack=True)
    
    #keep = np.where((ota_sky>0.0) & (ref_sky > 0.0) & (10000.0<ota_peak) & (10000.0<ref_peak) & (45000.0>ota_peak) & (45000.0>ref_peak) & (ota_fwhm < 100) & (ref_fwhm < 100))
    keep = np.where((ota_sky>0.0) & (ref_sky > 0.0) & (5000.0<ota_peak) & (5000.0<ref_peak) & (45000.0>ota_peak) & (45000.0>ref_peak) & (ota_fwhm < 100) & (ref_fwhm < 100))

    magA =  ota_mag[keep[0]]
    magRef =  ref_mag[keep[0]]
    
    print 'Calculating scaling factor for', img
    fitsref = odi.fits.open(ref_img)
    hduref = fitsref[0]
    fitsimg = odi.fits.open(img)
    hduimg = fitsimg[0]
    
    exp_ref = hduref.header['EXPTIME']
    exp_img = hduimg.header['EXPTIME']
    
    n = 1
    expRatio = (float(exp_ref)/float(exp_img))
    
    # use photometry over the whole image, so read in all the phot output files and append
    rat = np.power(10.0,-0.4*(magA-magRef))/expRatio
    if verbose:
        for i,r in enumerate(rat):
            print img, i, r
          
    sigTest = np.std(rat)
    print len(rat), np.mean(rat), np.median(rat), np.std(rat), n
    scale = np.mean(rat)
    std = np.std(rat)
    #if sigTest <= sigThreshold:
        #scale = np.mean(rat)
        #std = np.std(rat)
    #else:
        #while sigTest > sigThreshold:
            #magTempA = magA
            #magTempRef = magRef
            #magA = magTempA[np.where(abs(rat-np.median(rat))<sigTest)]
            #magRef = magTempRef[np.where(abs(rat-np.median(rat))<sigTest)]
            #rat = np.power(10.0,-0.4*(magA-magRef))/expRatio
            #for i,r in enumerate(rat):
                #print magA[i], magRef[i], r
            #sigTest = np.std(rat)
            
            #n = n + 1
            #if n > 10:
                #print "Iteration did not converge to sigma <", repr(sigThreshold),"for", img
                #print "Quitting..."
                #exit()
            #print len(rat), np.mean(rat), np.median(rat), np.std(rat), n
        #scale[img] = np.mean(rat)
        #std[img] = np.std(rat)
	#scale = np.mean(rat)
        #std = np.std(rat)
    dither_directory = odi.matchpath+'D'+img.split('.')[1][0]+filter
    if not os.path.exists(dither_directory):
	print 'Creating directory for dither',filter,'D'+img.split('.')[1][0]
	os.makedirs(dither_directory)
    dither_path = dither_directory+'/'
    print dither_path
    files_to_move = glob.glob('match/*')
    for file in files_to_move:
	if not file.startswith('match/D'):
	    shutil.move(file,dither_path+file.replace('match/',''))
    
    return scale, std


def main():
    pass

if __name__ == '__main__':
    main()
