import sys, os, glob, string
import numpy as np
import matplotlib.pyplot as plt
from pyraf import iraf
from tqdm import tqdm
import odi_config as odi
import pandas as pd
from astropy.wcs import WCS
from astropy.table import Table
from astropy.io import fits
from collections import OrderedDict

def get_sdss_coords(img, ota, inst,output='test.sdss'):
    formats = ['csv','xml','html']

    astro_url='http://skyserver.sdss3.org/public/en/tools/search/x_sql.aspx'
    public_url='http://skyserver.sdss.org/dr12/SkyserverWS/ImagingQuery/Cone?'

    default_url=public_url
    default_fmt='csv'

    hdulist = odi.fits.open(img.f)
    hdu = odi.tan_header_fix(hdulist[ota])
    
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

        #qry = "limit=5000&format=csv&imgparams=ra,dec,u,err_u,g,err_g,r,err_r,i,err_i,z,err_z,probPSF&specparams=none&ra="+repr(rac)+"&dec="+repr(decc)+"&radius="+repr(sizeam)+"&magType=psf"
        qry = "limit=5000&format=csv&imgparams=ra,dec,psfMag_u,psfMagErr_u,psfMag_g,psfMagErr_g,psfMag_r,psfMagErr_r,psfMag_i,psfMagErr_i,psfMag_z,psfMagErr_z,probPSF&specparams=none&ra="+repr(rac)+"&dec="+repr(decc)+"&radius="+repr(sizeam)+"&magType=psf"

        #print 'with query\n-->', qry
        print 'fetching SDSS sources around',rac,decc,'with radius',sizeam,'arcmin'
        url = default_url
        fmt = default_fmt
        writefirst = 1
        verbose = 0

        ofp = open(output,'w+')
        if verbose:
            odi.write_header(ofp,'#',url,qry)
        file_ = odi.httpquery(qry,url,fmt)
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

def refetch_sdss_coords(img, ota, gapmask, inst,gmaglim=19.,offline = False,source='sdss'):

    image = odi.reprojpath+'reproj_'+ota+'.'+img.stem()
    outcoords = odi.coordspath+'reproj_'+ota+'.'+img.base()+'.sdss'
    hdulist = odi.fits.open(image)
    
    hdu = odi.tan_header_fix(hdulist[0])
    xdim = hdu.header['NAXIS1']
    ydim = hdu.header['NAXIS2']

    if offline == False:
        formats = ['csv','xml','html']
        astro_url='http://skyserver.sdss3.org/public/en/tools/search/x_sql.aspx'
        public_url='http://skyserver.sdss.org/dr12/SkyserverWS/ImagingQuery/Cone?'
        default_url=public_url
        default_fmt='csv'
        if not os.path.isfile(outcoords):
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

            #qry = "limit=5000&format=csv&imgparams=ra,dec,u,err_u,g,err_g,r,err_r,i,err_i,z,err_z,probPSF&specparams=none&ra="+repr(rac)+"&dec="+repr(decc)+"&radius="+repr(sizeam)+"&magType=psf"
            qry = "limit=5000&format=csv&imgparams=ra,dec,psfMag_u,psfMagErr_u,psfMag_g,psfMagErr_g,psfMag_r,psfMagErr_r,psfMag_i,psfMagErr_i,psfMag_z,psfMagErr_z,probPSF&specparams=none&ra="+repr(rac)+"&dec="+repr(decc)+"&radius="+repr(sizeam)+"&magType=psf"


            #print 'with query\n-->', qry
            print 'fetching SDSS sources around',rac,decc,'with radius',sizeam,'arcmin'
            url = default_url
            fmt = default_fmt
            writefirst = 1
            verbose = 0

            ofp = open(outcoords,'w+')
            if verbose:
                odi.write_header(ofp,'#',url,qry)
            file_ = odi.httpquery(qry,url,fmt)
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
        w = odi.WCS(hdu.header)
        with open(odi.coordspath+'reproj_'+ota+'.'+img.base()+'.sdssxy', 'w+') as fxy:
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
                            print >> fxy, pixcrd2[0][0], pixcrd2[0][1], ras[i],decs[i],psfMag_u[i],psfMagErr_u[i],psfMag_g[i],psfMagErr_g[i],psfMag_r[i],psfMagErr_r[i],psfMag_i[i],psfMagErr_i[i],psfMag_z[i],psfMagErr_z[i]


    if offline == True:
        if source == 'sdss':
            outcoords =  odi.sdsspath+'offline_'+ota+'.'+img.base()+'.sdss'
            ras,decs,psfMag_u,psfMagErr_u,psfMag_g,psfMagErr_g,psfMag_r,psfMagErr_r,psfMag_i,psfMagErr_i,psfMag_z,psfMagErr_z = np.loadtxt(outcoords,usecols=(0,1,2,3,4,5,6,7,8,9,10,11), unpack=True, delimiter=',', skiprows=1)
        if source == 'twomass':
            outcoords =  odi.twomasspath+'offline_'+ota+'.'+img.base()+'.mass'
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
        if source == 'gaia':
            outcoords = odi.gaiapath+'offline_'+ota+'.'+img.base()+'.gaia'
            ras,decs = np.loadtxt(outcoords,usecols=(0,1), unpack=True, delimiter=',', skiprows=1)
            # Just creating dummy variables so that the file formats remain the same
            # for other functions
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
        tqdm.write('Using Ra and Dec from {:s} for reproject'.format(outcoords))
        w = odi.WCS(hdu.header)
        with open(odi.coordspath+'reproj_'+ota+'.'+img.base()+'.sdssxy', 'w+') as fxy:
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

def repoxy_offline(img, ota, gapmask, inst,gmaglim=19.,source='sdss'):
    image = odi.reprojpath+'reproj_'+ota+'.'+img.stem()
    hdulist = odi.fits.open(image)
    hdu = odi.tan_header_fix(hdulist[0])
    xdim = hdu.header['NAXIS1']
    ydim = hdu.header['NAXIS2']
    if source == 'sdss':
        outcoords =  odi.sdsspath+'offline_'+ota+'.'+img.base()+'.sdss'
        ras,decs,psfMag_u,psfMagErr_u,psfMag_g,psfMagErr_g,psfMag_r,psfMagErr_r,psfMag_i,psfMagErr_i,psfMag_z,psfMagErr_z = np.loadtxt(outcoords,usecols=(0,1,2,3,4,5,6,7,8,9,10,11), unpack=True, delimiter=',', skiprows=1)
        outputxy =  odi.coordspath+'reproj_'+ota+'.'+img.base()+'.sdssxy'
    if source == 'twomass':
        outcoords =  odi.twomasspath+'offline_'+ota+'.'+img.base()+'.mass'
        outputxy =  odi.coordspath+'reproj_'+ota+'.'+img.base()+'.massxy'
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
    if source == 'gaia':
        outcoords = odi.gaiapath+'offline_'+ota+'.'+img.base()+'.gaia'
        outputxy =  odi.coordspath+'reproj_'+ota+'.'+img.base()+'.gaiaxy'
        ras,decs = np.loadtxt(outcoords,usecols=(0,1), unpack=True, delimiter=',', skiprows=1)
        # Just creating dummy variables so that the file formats remain the same
        # for other functions
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
    tqdm.write('Using Ra and Dec from {:s} for reproject'.format(outcoords))
    w = odi.WCS(hdu.header)
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

def sdss_coords_full(img, inst,gmaglim=19.):
    formats = ['csv','xml','html']

    astro_url='http://skyserver.sdss3.org/public/en/tools/search/x_sql.aspx'
    #public_url='http://skyserver.sdss.org/SkyserverWS/dr12/ImagingQuery/Cone?'
    public_url='http://skyserver.sdss.org/dr12/SkyserverWS/ImagingQuery/Cone?'

    default_url=public_url
    default_fmt='csv'

    image = img
    outcoords = img.nofits()+'.sdss'

    hdulist = odi.fits.open(image)
    hdu = odi.tan_header_fix(hdulist[0])
    data = hdu.data
    # hdu = hdulist[ota]

    xdim = hdu.header['NAXIS1']
    ydim = hdu.header['NAXIS2']

    if not os.path.isfile(outcoords):
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

      #qry = "limit=5000&format=csv&imgparams=ra,dec,u,err_u,g,err_g,r,err_r,i,err_i,z,err_z,probPSF&specparams=none&ra="+repr(rac)+"&dec="+repr(decc)+"&radius="+repr(sizeam)+"&magType=psf"
      qry = "limit=10000&format=csv&imgparams=ra,dec,psfMag_u,psfMagErr_u,psfMag_g,psfMagErr_g,psfMag_r,psfMagErr_r,psfMag_i,psfMagErr_i,psfMag_z,psfMagErr_z,probPSF&specparams=none&ra="+repr(rac)+"&dec="+repr(decc)+"&radius="+repr(sizeam)+"&magType=psf"

      #print 'with query\n-->', qry
      print 'fetching SDSS sources around',rac,decc,'with radius',sizeam,'arcmin'
      url = default_url
      fmt = default_fmt
      writefirst = 1
      verbose = 0

      ofp = open(outcoords,'w+')
      if verbose:
          odi.write_header(ofp,'#',url,qry)
      file_ = odi.httpquery(qry,url,fmt)
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

    w = odi.WCS(hdu.header)
    with open(img.nofits()+'.wcs.coo','w+') as f:
        with open(img.nofits()+'.sdssxy', 'w+') as fxy:
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

# def get_sdss_coords_offline(img, ota, inst,output='test.sdss'):
#     hdulist = odi.fits.open(img.f)
# 
#     hdu = odi.tan_header_fix(hdulist[ota])
#     xdim = hdu.header['NAXIS1']
#     ydim = hdu.header['NAXIS2']
# 
#     sdss_cat_img = hdulist['CAT.PHOTCALIB']
#     sdss_cat_img_df = pd.DataFrame.from_dict(sdss_cat_img.data)
#     hdulist.close()
# 
#     ota = float(ota.strip('OTA.SCI'))
#     ota_matches_df = sdss_cat_img_df.iloc[np.where(sdss_cat_img_df['ODI_OTA'] == ota)]
# 
#     needed_columns = ['SDSS_RA','SDSS_DEC','SDSS_MAG_U',
#                       'SDSS_ERR_U', u'SDSS_MAG_G', u'SDSS_ERR_G', u'SDSS_MAG_R',
#                       'SDSS_ERR_R', u'SDSS_MAG_I', u'SDSS_ERR_I', u'SDSS_MAG_Z',
#                       'SDSS_ERR_Z','ODI_OTA','ODI_OTA']
# 
#     output_df = ota_matches_df[needed_columns]
#     output_df.to_csv(output,index=False)
#     return xdim, ydim

def refetch_sdss_coords_offline(img, ota, gapmask, inst,gmaglim=19.):
    image = odi.reprojpath+'reproj_'+ota+'.'+img.stem()
    outcoords = odi.coordspath+'reproj_'+ota+'.'+img.base()+'.sdss'

    hdulist = odi.fits.open(image)
    hdu = odi.tan_header_fix(hdulist[0])
    xdim = hdu.header['NAXIS1']
    ydim = hdu.header['NAXIS2']

    if not os.path.isfile(outcoords):
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



    ras,decs,psfMag_u,psfMagErr_u,psfMag_g,psfMagErr_g,psfMag_r,psfMagErr_r,psfMag_i,psfMagErr_i,psfMag_z,psfMagErr_z = np.loadtxt(outcoords,usecols=(0,1,2,3,4,5,6,7,8,9,10,11), unpack=True, delimiter=',', skiprows=2)
    probPSF = np.loadtxt(outcoords, usecols=(12,), dtype=int, unpack=True, delimiter=',', skiprows=2)
    # print ras, decs

    w = odi.WCS(hdu.header)

    with open(odi.coordspath+'reproj_'+ota+'.'+img.base()+'.sdssxy', 'w+') as fxy:
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

def get_sdss_coords_offline(img, ota, inst,output='test.sdss'):
    """
    Pull out and parse the ``CAT.PHOTCALIB`` table from ``img`` header. This
    function will separate the SDSS stars in ``CAT.PHOTCALIB`` based on which
    ``ota`` they fall on.

    Parameters
    ----------
    img : str
        Name of image
    ota : str
        Name of OTA
    int : str
        Version of ODI used, ``podi`` or ``5odi``

    Returns
    -------
    xdim : int
        Size of OTA in the x direction ``NAXIS1``
    ydim : int
        Size of OTA in the y direction ``NAXIS2``

    Note
    ----
    If the images being processed do not fall in the SDSS footprint,
    the QuickReduce pipeline will use PanSTARRS. This function will still
    pull out these stars and treat them as SDSS stars. There will be no
    ``u`` magnitudes available, however.
    """
    hdulist = odi.fits.open(img.f)
    hdu = odi.tan_header_fix(hdulist[ota])

    xdim = hdu.header['NAXIS1']
    ydim = hdu.header['NAXIS2']
    try:
        sdss_cat_img = hdulist['CAT.PHOTCALIB']
        cat_img_data = Table.read(sdss_cat_img, format='fits')
        # print cat_img_data.colnames
        # force little-endian byte order to make FITS play nice with pandas
        sdss_cat_img_df = cat_img_data.to_pandas()
        # sdss_cat_img_df = pd.DataFrame.from_dict(cat_img_dict)
        # print sdss_cat_img_df.keys()
        ota = float(ota.strip('OTA.SCI'))
        # print 'catalog source:', hdulist[0].header['PHOTMCAT']
        if 'sdss_dr' in hdulist[0].header['PHOTMCAT']:
            try:
                # print sdss_cat_img_df.columns
                ota_matches_df = sdss_cat_img_df.iloc[np.where(sdss_cat_img_df[u'ODI_OTA'] == ota)]
                needed_columns = [u'REF_RA',u'REF_DEC',u'REF_U',
                                  u'REF_ERR_U', u'REF_G', u'REF_ERR_G', u'REF_R',
                                  u'REF_ERR_R', u'REF_I', u'REF_ERR_I', u'REF_Z',
                                  u'REF_ERR_Z', u'ODI_OTA']
    
                output_df = ota_matches_df[needed_columns]
                output_df.to_csv(output,index=False)
            except KeyError:
                oditable = hdulist['CAT.ODI'].data
                oditalbe_df = pd.DataFrame.from_dict(oditable)
    
                ODI_RA = np.squeeze(np.array(oditalbe_df['RA']))
                ODI_DEC = np.squeeze( np.array(oditalbe_df['DEC']))
                ODI_OTA = np.squeeze( np.array(oditalbe_df['OTA']))
    
                junkdict = OrderedDict([(u'ODI_RA',ODI_RA),
                                        (u'ODI_DEC',ODI_DEC),
                                        (u'ODI_OTA',ODI_OTA.astype(float))])
                junk_df = pd.DataFrame.from_dict(junkdict)
    
                matched_df = pd.merge(sdss_cat_img_df,junk_df ,on = [u'ODI_RA',u'ODI_DEC'],how='inner')
                # print matched_df.columns
                needed_columns = np.insert(sdss_cat_img_df.columns.values,0,u'ODI_OTA')
    
                full_df = matched_df[needed_columns]
                ota_matches_df = full_df.iloc[np.where(full_df[u'ODI_OTA'] == ota)]
                needed_columns = [u'SDSS_RA',u'SDSS_DEC',
                                  u'SDSS_MAG_U',u'SDSS_ERR_U',
                                  u'SDSS_MAG_G', u'SDSS_ERR_G',
                                  u'SDSS_MAG_R',u'SDSS_ERR_R',
                                  u'SDSS_MAG_I', u'SDSS_ERR_I',
                                  u'SDSS_MAG_Z',u'SDSS_ERR_Z',
                                  u'ODI_OTA']
                output_df = ota_matches_df[needed_columns]
                output_df.to_csv(output,index=False)
        else:
            ota_matches_df = sdss_cat_img_df.iloc[np.where(sdss_cat_img_df[u'ODI_OTA'] == ota)]
            ota_matches_df = ota_matches_df.reset_index()
            junk_u = np.ones(len(ota_matches_df))
            junk_u_err = np.ones(len(ota_matches_df))
            ota_matches_df[u'IPP_MAG_U'] = junk_u
            ota_matches_df[u'IPP_ERR_U'] = junk_u_err
    
            needed_columns = [u'IPP_RA', u'IPP_DEC',
                              u'IPP_MAG_U', u'IPP_ERR_U',
                              u'IPP_MAG_G', u'IPP_ERR_G',
                              u'IPP_MAG_R', u'IPP_ERR_R',
                              u'IPP_MAG_I', u'IPP_ERR_I',
                              u'IPP_MAG_Z',u'IPP_ERR_Z',
                              u'ODI_OTA']
    
            output_df = ota_matches_df[needed_columns]
            output_df.to_csv(output,index=False)
    
        if 'SDSS' in hdulist[0].header['PHOTMCAT']:
            try:
                # print sdss_cat_img_df.columns
                ota_matches_df = sdss_cat_img_df.iloc[np.where(sdss_cat_img_df[u'ODI_OTA'] == ota)]
                needed_columns = [u'SDSS_RA',u'SDSS_DEC',u'SDSS_MAG_U',
                                  u'SDSS_ERR_U', u'SDSS_MAG_G', u'SDSS_ERR_G', u'SDSS_MAG_R',
                                  u'SDSS_ERR_R', u'SDSS_MAG_I', u'SDSS_ERR_I', u'SDSS_MAG_Z',
                                  u'SDSS_ERR_Z', u'ODI_OTA']
    
                output_df = ota_matches_df[needed_columns]
                output_df.to_csv(output,index=False)
            except KeyError:
                oditable = hdulist['CAT.ODI'].data
                oditalbe_df = pd.DataFrame.from_dict(oditable)
    
                ODI_RA = np.squeeze(np.array(oditalbe_df['RA']))
                ODI_DEC = np.squeeze( np.array(oditalbe_df['DEC']))
                ODI_OTA = np.squeeze( np.array(oditalbe_df['OTA']))
    
                junkdict = OrderedDict([(u'ODI_RA',ODI_RA),
                                        (u'ODI_DEC',ODI_DEC),
                                        (u'ODI_OTA',ODI_OTA.astype(float))])
                junk_df = pd.DataFrame.from_dict(junkdict)
    
                matched_df = pd.merge(sdss_cat_img_df,junk_df ,on = [u'ODI_RA',u'ODI_DEC'],how='inner')
                # print matched_df.columns
                needed_columns = np.insert(sdss_cat_img_df.columns.values,0,u'ODI_OTA')
    
                full_df = matched_df[needed_columns]
                ota_matches_df = full_df.iloc[np.where(full_df[u'ODI_OTA'] == ota)]
                needed_columns = [u'SDSS_RA',u'SDSS_DEC',
                                  u'SDSS_MAG_U',u'SDSS_ERR_U',
                                  u'SDSS_MAG_G', u'SDSS_ERR_G',
                                  u'SDSS_MAG_R',u'SDSS_ERR_R',
                                  u'SDSS_MAG_I', u'SDSS_ERR_I',
                                  u'SDSS_MAG_Z',u'SDSS_ERR_Z',
                                  u'ODI_OTA']
                output_df = ota_matches_df[needed_columns]
                output_df.to_csv(output,index=False)
    except KeyError:
        print 'missing PHOTCALIB table, skipping SDSS'
    hdulist.close()
    return xdim, ydim

def get_gaia_coords(img,ota,inst,output='test.gaia',cluster=False,**kwargs):
    """
    Query the online Gaia DR1 based on the central coordinates of the current
    OTA. If the ``cluster`` flag is set to ``True``, the querey will avoid
    a crowded region based on coordinates and a radius set by the user in
    the configuration files.

    Parameters
    ----------
    img : ODIImage or StackedImage object
        Name of image
    ota : str
        Name of OTA
    int : str
        Version of ODI used, ``podi`` or ``5odi``

    """
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    try:
        from astroquery.vizier import Vizier
        from astropy import __version__ as astropyversion
    except ImportError:
        print "astroquery not installed"
        print "try  pip --user --no-deps install astroquery or contact admin"
    hdulist = fits.open(img.f)
    if ota=='None':
        hdu_ota = hdulist[0]

    else:
        hdu_ota = odi.tan_header_fix(hdulist[ota])

    w = WCS(hdu_ota.header)

    naxis1 = hdu_ota.header['NAXIS1']
    naxis2 = hdu_ota.header['NAXIS2']
    ota_center_radec = w.wcs_pix2world([[naxis1/2.,naxis2/2.]],1)

    corners = w.calc_footprint()

    center_skycoord = SkyCoord(ota_center_radec[0][0]*u.deg,
                               ota_center_radec[0][1]*u.deg,frame='icrs')
    corner_skycoord = SkyCoord(corners[0,0]*u.deg,
                               corners[0,1]*u.deg,frame='icrs')
    cone_radius = center_skycoord.separation(corner_skycoord).value
    # tqdm.write('{:4.0f} {:4.0f} {:6.4f}'.format(naxis1/2., naxis2/2., cone_radius))

    #Set up vizier query for Gaia DR1
    #Taken from example at: github.com/mommermi/photometrypipeline
    vquery = Vizier(columns=['RA_ICRS', 'DE_ICRS',
                             'e_RA_ICRS', 'e_DE_ICRS',
                             'phot_g_mean_mag'],
                    column_filters={"phot_g_mean_mag":
                                    ("<%f" % 21.0)},
                    row_limit = -1)
    gaia_table = vquery.query_region(SkyCoord(ra=ota_center_radec[0][0],
                                              dec=ota_center_radec[0][1],
                                              unit=(u.deg, u.deg),
                                              frame='icrs'),
                                     radius=cone_radius*u.deg,
                                     catalog=['I/337/gaia'])[0]

    # print gaia_table
    
    hdulist.close()
    if cluster == True:
        try:
            racenter = kwargs['racenter']
            deccenter = kwargs['deccenter']
            min_radius = kwargs['min_radius']
            G_lim = kwargs['G_lim']
        except KeyError:
            print 'Must provide racenter, deccenter, and min_radius'
        cluster_center = SkyCoord(racenter*u.degree
                                  ,deccenter*u.degree,
                                  frame='icrs')
        gaia_coords = SkyCoord(gaia_table['RA_ICRS'],
                               gaia_table['DE_ICRS'],frame='icrs')
        dist_from_center = cluster_center.separation(gaia_coords).arcmin
        gaia_table['dis'] = dist_from_center
        # ota_gaia_df = ota_gaia_df[(ota_gaia_df.dis >= min_radius) &
                                #   (ota_gaia_df.phot_g_mean_mag <= G_lim)]
        gaia_table = gaia_table[gaia_table['dis'] > min_radius]

    ra_min, ra_max = min(corners[:,0]), max(corners[:,0])
    dec_min, dec_max = min(corners[:,1]), max(corners[:,1])
    # print ra_min, ra_max, dec_min, dec_max

    gaia_table_cut = gaia_table[(gaia_table['RA_ICRS'] > ra_min) &
                                (gaia_table['RA_ICRS'] < ra_max) &
                                (gaia_table['DE_ICRS'] > dec_min) &
                                (gaia_table['DE_ICRS'] < dec_max)]
    gaia_table_cut['e_RA_ICRS'].convert_unit_to(u.deg)
    gaia_table_cut['e_DE_ICRS'].convert_unit_to(u.deg)

    ota_gaia_df = gaia_table_cut.to_pandas()

    cols_needed = ['RA_ICRS','DE_ICRS','__Gmag_','e_RA_ICRS','e_DE_ICRS']

    ota_gaia_df = ota_gaia_df[cols_needed]
    ota_gaia_df.columns = ['ra', 'dec','phot_g_mean_mag','e_ra','e_dec']

    gaia_catalog_out = output
    ota_gaia_df.to_csv(gaia_catalog_out,
                       columns=['ra', 'dec','phot_g_mean_mag','e_ra','e_dec'],
                       index=False)
    return ota_gaia_df


if __name__ == '__main__':
    from odi_config import ODIImage
    img = ODIImage("20161029T215338.1_AGC748738_odi_g.7498.fits", 1, '5odi')
    get_gaia_coords(img, img.otas.values()[0], img.inst, output='test.gaia', cluster=False)
