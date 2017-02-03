from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.wcs import WCS
#from rand_bkg import bkg_boxes
from astropy.convolution import Gaussian2DKernel
from astropy.stats import sigma_clipped_stats
from photutils import detect_sources
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
from astropy.table import Table

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

    sdss_cat_img = hdulist['CAT.PHOTCALIB']
    cat_img_data = Table.read(sdss_cat_img, format='fits')
    # print cat_img_data.colnames
    # force little-endian byte order to make FITS play nice with pandas
    sdss_cat_img_df = cat_img_data.to_pandas()
    # sdss_cat_img_df = pd.DataFrame.from_dict(cat_img_dict)
    # print sdss_cat_img_df.keys()
    ota = float(ota.strip('OTA.SCI'))
    print 'catalog source:', hdulist[0].header['PHOTMCAT']
    if 'sdss' in hdulist[0].header['PHOTMCAT'] or 'SDSS' in hdulist[0].header['PHOTMCAT']:
        try:
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

            junkdict = OrderedDict([('ODI_RA',ODI_RA),
                                    ('ODI_DEC',ODI_DEC),
                                    ('ODI_OTA',ODI_OTA.astype(float))])
            junk_df = pd.DataFrame.from_dict(junkdict)

            matched_df = pd.merge(sdss_cat_img_df,junk_df ,on = ['ODI_RA','ODI_DEC'],how='inner')

            needed_columns = np.insert(sdss_cat_img_df.columns.values,0,'ODI_OTA')

            full_df = matched_df[needed_columns]
            ota_matches_df = full_df.iloc[np.where(full_df['ODI_OTA'] == ota)]
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
        ota_matches_df = sdss_cat_img_df.iloc[np.where(sdss_cat_img_df['ODI_OTA'] == ota)]
        ota_matches_df = ota_matches_df.reset_index()
        junk_u = np.ones(len(ota_matches_df))
        junk_u_err = np.ones(len(ota_matches_df))
        ota_matches_df['IPP_MAG_U'] = junk_u
        ota_matches_df['IPP_ERR_U'] = junk_u_err

        needed_columns = ['IPP_RA', 'IPP_DEC',
                          'IPP_MAG_U', 'IPP_ERR_U',
                          'IPP_MAG_G', 'IPP_ERR_G',
                          'IPP_MAG_R', 'IPP_ERR_R',
                          'IPP_MAG_I', 'IPP_ERR_I',
                          'IPP_MAG_Z','IPP_ERR_Z',
                          'ODI_OTA']

        output_df = ota_matches_df[needed_columns]
        output_df.to_csv(output,index=False)

    hdulist.close()
    return xdim, ydim

def get_2mass_coords_offline(img, ota, inst,output='test.mass'):
    """
    Pull out and parse the ``CAT.ODI+2MASS`` table from ``img`` header. This
    function will separate the 2MASS stars in ``CAT.ODI+2MASS`` based on which
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
    """
    hdulist = odi.fits.open(img.f)
    hdu = odi.tan_header_fix(hdulist[ota])
    
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
    """
    Query the online Gaia DR1 based on the central coordinates of the current
    OTA. If the ``cluster`` flag is set to ``True``, the querey will avoid
    a crowded region based on coordinates and a radius set by the user in
    the configuration files.

    Parameters
    ----------
    img : str
        Name of image
    ota : str
        Name of OTA
    int : str
        Version of ODI used, ``podi`` or ``5odi``

    """
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    hdulist = fits.open(img.f)
    hdu_ota = odi.tan_header_fix(hdulist[ota])
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
            G_lim = kwargs['G_lim']
        except KeyError:
            print 'Must provide racenter, deccenter, and min_radius'
        cluster_center = SkyCoord(racenter*u.degree
                                  ,deccenter*u.degree,
                                  frame='icrs')
        gaia_coords = SkyCoord(ota_gaia_df.ra.values*u.deg,
                               ota_gaia_df.dec.values*u.deg,frame='icrs')
        dist_from_center = cluster_center.separation(gaia_coords).arcmin
        ota_gaia_df['dis'] = dist_from_center
        ota_gaia_df = ota_gaia_df[(ota_gaia_df.dis >= min_radius) &
                                  (ota_gaia_df.phot_g_mean_mag <= G_lim)]

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
