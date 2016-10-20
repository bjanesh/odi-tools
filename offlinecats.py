from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.wcs import WCS
#from rand_bkg import bkg_boxes
from astropy.convolution import Gaussian2DKernel
from astropy.stats import sigma_clipped_stats
from photutils.detection import detect_sources
from photutils.utils import random_cmap
from scipy.ndimage import binary_dilation
import pandas as pd
import numpy as np
import os
import glob
from collections import OrderedDict
from astropy.io import fits
import odi_config as odi
from gaia.tap import cone_search

def get_sdss_coords_offline(img, ota, inst,output='test.sdss'):
    hdulist = odi.fits.open(img)
    hdu = hdulist[ota]

    if inst == 'podi':
        pvlist = hdu.header['PV*']
        for pv in pvlist:
            tpv = 'T'+pv
            hdu.header.rename_keyword(pv, tpv, force=False)
    xdim = hdu.header['NAXIS1']
    ydim = hdu.header['NAXIS2']

    sdss_cat_img = hdulist['CAT.PHOTCALIB']
    sdss_cat_img_df = pd.DataFrame.from_dict(sdss_cat_img.data)

    ota = float(ota.strip('OTA.SCI'))
    try:
        ota_matches_df = sdss_cat_img_df.iloc[np.where(sdss_cat_img_df['ODI_OTA'] == ota)]
        needed_columns = ['SDSS_RA','SDSS_DEC','SDSS_MAG_U',
                          'SDSS_ERR_U', u'SDSS_MAG_G', u'SDSS_ERR_G', u'SDSS_MAG_R',
                          'SDSS_ERR_R', u'SDSS_MAG_I', u'SDSS_ERR_I', u'SDSS_MAG_Z',
                          'SDSS_ERR_Z','ODI_OTA']

        output_df = ota_matches_df[needed_columns]
        output_df.to_csv(output,index=False)
    except KeyError:
        oditable = hdulist['CAT.ODI'].data
        oditalbe_df = pd.DataFrame.from_dict(oditable)

        ODI_RA = np.squeeze(np.array(oditalbe_df['RA']))
        ODI_DEC = np.squeeze( np.array(oditalbe_df['DEC']))
        ODI_OTA = np.squeeze( np.array(oditalbe_df['OTA']))

        junkdict = OrderedDict([('ODI_RA',ODI_RA),
                                ('ODI_DEC',ODI_DEC),
                                ('ODI_OTA',ODI_OTA.astype(float))])
        junk_df = pd.DataFrame.from_dict(junkdict)

        matched_df = pd.merge(sdss_cat_img_df,junk_df ,on = ['ODI_RA','ODI_DEC'],how='inner')

        needed_columns = np.insert(sdss_cat_img_df.columns.values,0,'ODI_OTA')

        full_df = matched_df[needed_columns]
        ota_matches_df = full_df.iloc[np.where(full_df['ODI_OTA'] == ota)]
        needed_columns = ['SDSS_RA','SDSS_DEC','SDSS_MAG_U',
                          'SDSS_ERR_U', u'SDSS_MAG_G', u'SDSS_ERR_G', u'SDSS_MAG_R',
                          'SDSS_ERR_R', u'SDSS_MAG_I', u'SDSS_ERR_I', u'SDSS_MAG_Z',
                          'SDSS_ERR_Z','ODI_OTA']
        output_df = ota_matches_df[needed_columns]
        output_df.to_csv(output,index=False)
    hdulist.close()
    return xdim, ydim

def get_2mass_coords_offline(img, ota, inst,output='test.mass'):
    hdulist = odi.fits.open(img)
    hdu = hdulist[ota]

    if inst == 'podi':
        pvlist = hdu.header['PV*']
        for pv in pvlist:
            tpv = 'T'+pv
            hdu.header.rename_keyword(pv, tpv, force=False)
    xdim = hdu.header['NAXIS1']
    ydim = hdu.header['NAXIS2']

    mass_cat_img = hdulist['CAT.ODI+2MASS']
    mass_cat_img_df = mass_cat_img.data
    hdulist.close()

    ota = float(ota.strip('OTA.SCI'))
    ota_matches_df = mass_cat_img_df[np.where(mass_cat_img_df['OTA'] == ota)]

    junk_df = pd.DataFrame.from_dict(ota_matches_df)

    output_df = junk_df
    output_df.to_csv(output,index=False)
    return xdim, ydim

def get_gaia_coords(img,ota,inst,output='test.gaia',cluster=False,**kwargs):
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    hdulist = fits.open(img)
    hdu_ota = hdulist[ota]
    w = WCS(hdu_ota.header)
    ota_center_radec = w.wcs_pix2world([[2018.0,2007.5]],1)

    coord1 = w.wcs_pix2world([[1,1]],1)
    coord2 = w.wcs_pix2world([[4036,1]],1)
    coord3 = w.wcs_pix2world([[4036,4015]],1)
    coord4 = w.wcs_pix2world([[1,4015]],1)
    center_skycoord = SkyCoord(ota_center_radec[0][0]*u.deg,
                               ota_center_radec[0][1]*u.deg,frame='icrs')
    corner_skycoord = SkyCoord(coord2[0][0]*u.deg,
                               coord2[0][1]*u.deg,frame='icrs')
    cone_radius = center_skycoord.separation(corner_skycoord).value

    print 'Retrieving Gaia sources for: ', ota
    ota_gaia_sources = cone_search(ota_center_radec[0][0],
                                   ota_center_radec[0][1],
                                   cone_radius)
    ota_gaia_df = ota_gaia_sources.to_pandas()
    hdulist.close()
    if cluster == True:
        try:
            racenter = kwargs['racenter']
            deccenter = kwargs['deccenter']
            min_radius = kwargs['min_radius']
        except KeyError:
            print 'Must provide racenter, deccenter, and min_radius'
        cluster_center = SkyCoord(racenter*u.degree
                                  ,deccenter*u.degree,
                                  frame='icrs')
        gaia_coords = SkyCoord(ota_gaia_df.ra.values*u.deg,
                               ota_gaia_df.dec.values*u.deg,frame='icrs')
        dist_from_center = cluster_center.separation(gaia_coords).arcmin
        ota_gaia_df['dis'] = dist_from_center
        ota_gaia_df = ota_gaia_df[ota_gaia_df.dis >= min_radius]

    ra_min, ra_max = coord1[0][0],coord2[0][0]
    dec_min, dec_max = coord1[0][1],coord3[0][1]

    ota_gaia_df = ota_gaia_df[(ota_gaia_df.ra > ra_min) &
                              (ota_gaia_df.ra < ra_max) &
                              (ota_gaia_df.dec > dec_min) &
                              (ota_gaia_df.dec < dec_max)]

    gaia_catalog_out = output
    ota_gaia_df.to_csv(gaia_catalog_out,
                       columns=['ra', 'dec','phot_g_mean_mag'],
                       index=False)


def main():
    pass

if __name__ == '__main__':
    main()
