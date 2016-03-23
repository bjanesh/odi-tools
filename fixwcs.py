import sys, os, glob, string
import numpy as np
import astropy as ast
import matplotlib.pyplot as plt
from pyraf import iraf
from tqdm import tqdm
import odi_config as odi

def list_wcs_coords(img, ota, gapmask, inst,output='radec.coo', gmaglim=20., stars_only=True, offline = False, source = 'sdss'):
    # outputs a list of coordinates in the proper formatting to be used as input to 
    # iraf.msccmatch. RETURNS a list of pixel coordinates for plotting
    if offline == False:
        xdim, ydim = odi.get_sdss_coords(img, ota, inst,output=odi.coordspath+img[0:-5]+'.'+ota+'.sdss')
        ras,decs,psfMag_u,psfMagErr_u,psfMag_g,psfMagErr_g,psfMag_r,psfMagErr_r,psfMag_i,psfMagErr_i,psfMag_z,psfMagErr_z = np.loadtxt(odi.coordspath+img[0:-5]+'.'+ota+'.sdss',usecols=(0,1,2,3,4,5,6,7,8,9,10,11), unpack=True, delimiter=',', skiprows=2)
        probPSF = np.loadtxt(odi.coordspath+img[0:-5]+'.'+ota+'.sdss', usecols=(12,), dtype=int, unpack=True, delimiter=',', skiprows=2)
        coords2 = zip(ras[np.where((psfMag_g<gmaglim) & (probPSF==1))],decs[np.where((psfMag_g<gmaglim) & (probPSF==1))])
    if offline == True and source == 'sdss':
        sdss_cat = odi.sdsspath+'offline_'+ota+'.'+str(img[16:-5])+'.sdss'
        print 'Using Ra and Dec from:', sdss_cat,'for fixwcs'
        ras,decs,psfMag_u,psfMagErr_u,psfMag_g,psfMagErr_g,psfMag_r,psfMagErr_r,psfMag_i,psfMagErr_i,psfMag_z,psfMagErr_z = np.loadtxt(sdss_cat,usecols=(0,1,2,3,4,5,6,7,8,9,10,11), unpack=True, delimiter=',', skiprows=1)
        coords2 = zip(ras[np.where(psfMag_g<gmaglim)],decs[np.where(psfMag_g<gmaglim)])
    if offline == True and source == 'twomass':
        twomass_cat = odi.twomasspath+'offline_'+ota+'.'+str(img[16:-5])+'.mass'
        ras,decs = np.loadtxt(twomass_cat,usecols=(2,3), unpack=True, delimiter=',', skiprows=1)
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
        coords2 = zip(ras,decs)

    hdulist = odi.fits.open(img)
    hdu = hdulist[ota]
    if inst == 'podi':
        pvlist = hdu.header['PV*']
        for pv in pvlist:
            tpv = 'T'+pv
            hdu.header.rename_keyword(pv, tpv, force=False)
    if offline == True:
        xdim = hdu.header['NAXIS1']
        ydim = hdu.header['NAXIS2']
    
    w = odi.WCS(hdu.header)
    pixcrd2 = w.wcs_world2pix(coords2, 1)
    pixid = []
    with open(odi.coordspath+output,'w+') as f:
        with open(odi.coordspath+output+'.pix', 'w+') as fp:
            with open(odi.coordspath+img[0:-5]+'.'+ota+'.sdssxy', 'w+') as fxy:
                for i,c in enumerate(coords2):
                    if  20.0 <= pixcrd2[i,0] < xdim-100.0 and 20.0 <= pixcrd2[i,1] < ydim-100.0:
                        # make an image cutout of the gap mask
                        x, y = int(round(pixcrd2[i,0])), int(round(pixcrd2[i,1]))
                        cutout = gapmask[x-30:x+30,y-30:y+30]
                        # print cutout
                        if not (cutout.astype(bool)).any():
                            pixid.append(i)
                            r, d = odi.deg_to_sex(c[0], c[1])
                            print >> f, r, d, psfMag_g[i]
                            print >> fp, pixcrd2[i,0], pixcrd2[i,1], i, 'm'
                            print >> fxy, pixcrd2[i,0], pixcrd2[i,1], ras[i],decs[i],psfMag_u[i],psfMagErr_u[i],psfMag_g[i],psfMagErr_g[i],psfMag_r[i],psfMagErr_r[i],psfMag_i[i],psfMagErr_i[i],psfMag_z[i],psfMagErr_z[i]

    pixid = np.array(pixid)
    pixcrd3 = pixcrd2[pixid]
    hdulist.close()
    return pixcrd3

def fix_wcs(img, ota, coords='radec.coo', iters=3):
  image = odi.illcorpath+'illcor_'+ota+'.'+str(img[16:])
  iraf.mscred(_doprint=0)
  iraf.unlearn(iraf.mscred.msccmatch)
  # otaext = {'33':'[1]','34':'[2]','44':'[3]','43':'[4]','42':'[5]','32':'[6]','22':'[7]','23':'[8]','24':'[9]'}
  for i in range(iters):
      fix = iraf.msccmatch(input=image,coords=odi.coordspath+coords,usebpm='no',verbose='yes',nsearch=1000,search=10,rsearch=5.5,cfrac=.9,csig=1.0,nfit=3,rms=15,maxshif=50,fitgeom="general",update='yes',interac='no',fit='no',accept='yes', Stdout=1)
      print 'fixing WCS for',img+'['+ota+'], iter ='+repr(i)
      print fix[-6]
      print fix[-5]
      print fix[-4]
      print fix[-3]
      print fix[-2]
      
def fix_wcs_full(img, coords='radec.coo', iters=3):
    print coords 
    iraf.mscred(_doprint=0)
    iraf.unlearn(iraf.mscred.msccmatch)
    # otaext = {'33':'[1]','34':'[2]','44':'[3]','43':'[4]','42':'[5]','32':'[6]','22':'[7]','23':'[8]','24':'[9]'}
    for i in range(iters):
        fix = iraf.msccmatch(input=img,coords=coords,usebpm='no',verbose='yes',nsearch=1000,search=10,rsearch=5.5,cfrac=.9,csig=1.0,nfit=3,rms=15,maxshif=50,fitgeom="general",update='yes',interac='no',fit='no',accept='yes', Stdout=1)
        print 'fixing WCS for',img+', iter ='+repr(i)
        print fix[-6]
        print fix[-5]
        print fix[-4]
        print fix[-3]
        print fix[-2]	  

def main():
    pass

if __name__ == '__main__':
    main()
    