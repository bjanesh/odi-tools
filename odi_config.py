from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.io import fits
from astropy.wcs import WCS
from astropy.convolution import Gaussian2DKernel
from astropy.stats import sigma_clipped_stats
from photutils.segmentation import detect_sources
# from photutils.utils import random_cmap
from scipy.ndimage import binary_dilation

from odi_calibrate import *
from odi_illcor import *
from odi_helpers import *
from odi_coords import *
from odi_scale import *
from full_calibrate import *
from full_phot import *

podi_dictionary = {
    1: 'OTA33.SCI',
    2: 'OTA34.SCI',
    3: 'OTA44.SCI',
    4: 'OTA43.SCI',
    5: 'OTA42.SCI',
    6: 'OTA32.SCI',
    7: 'OTA22.SCI',
    8: 'OTA23.SCI',
    9: 'OTA24.SCI'
}

odi5_dictionary = {
    1 :'OTA33.SCI',
    2 :'OTA34.SCI',
    3 :'OTA43.SCI',
    4 :'OTA44.SCI',
    5 :'OTA32.SCI',
    6 :'OTA23.SCI',
    7 :'OTA24.SCI',
    8 :'OTA42.SCI',
    9 :'OTA35.SCI',
    10:'OTA53.SCI',
    11:'OTA45.SCI',
    12:'OTA54.SCI',
    13:'OTA22.SCI',
    14:'OTA25.SCI',
    15:'OTA52.SCI',
    16:'OTA55.SCI',
    17:'OTA31.SCI',
    18:'OTA13.SCI',
    19:'OTA41.SCI',
    20:'OTA14.SCI',
    21:'OTA36.SCI',
    22:'OTA46.SCI',
    23:'OTA21.SCI',
    24:'OTA12.SCI',
    25:'OTA15.SCI',
    26:'OTA51.SCI',
    27:'OTA26.SCI',
    28:'OTA56.SCI',
    29:'OTA11.SCI',
    30:'OTA16.SCI'
}

odi5mosaic_dictionary = {
    1 :'OTA33.SCI',
    2 :'OTA34.SCI',
    3 :'OTA43.SCI',
    4 :'OTA44.SCI',
    5 :'OTA32.SCI',
    6 :'OTA23.SCI',
    7 :'OTA24.SCI',
    8 :'OTA42.SCI',
    13:'OTA22.SCI',
}

# for simple demonstrations and tests, only use the central 4 OTAs for speed purposes
# valid for any configuration (podi, 5odi, mosaic)
test_dictionary = {
    1 :'OTA33.SCI',
    2 :'OTA34.SCI',
    3 :'OTA43.SCI',
    4 :'OTA44.SCI'
}

class ODIImage:
    def __init__(self, filename, dither, inst):
        self.f = filename
        self.d = dither
        self.inst = inst
        if self.inst == '5odi':
            self.otas = odi5_dictionary
        elif self.inst == 'podi':
            self.otas = podi_dictionary
        elif self.inst == 'mosaic':
            self.otas = odi5mosaic_dictionary
        else:
            raise ValueError('Instrument not recognized!')
        
    def nofits(self):
        return str(self.f[:-5])
        
    def dither(self):
        return repr(self.d)
        
    def stem(self):
        return repr(self.d)+str(self.f[17:])
        
    def base(self):
        return repr(self.d)+str(self.f[17:-5])
        
class StackedImage:
    def __init__(self, filename):
        self.f = filename

bpmdirectory = 'bpmasks'
if not os.path.exists(bpmdirectory):
    print('Creating directory for bad pixel masks...')
    os.makedirs(bpmdirectory)

bppath = bpmdirectory+'/'

illcordirectory = 'illcor'
if not os.path.exists(illcordirectory):
    print('Creating directory for illumination corrected ota images...')
    os.makedirs(illcordirectory)

illcorpath = illcordirectory+'/'

reprojdirectory = 'reproj'
if not os.path.exists(reprojdirectory):
    print('Creating directory for reprojected ota images...')
    os.makedirs(reprojdirectory)

reprojpath = reprojdirectory+'/'

bgsubdirectory = 'bgsub'
if not os.path.exists(bgsubdirectory):
    print('Creating directory for background subtracted ota images...')
    os.makedirs(bgsubdirectory)

bgsubpath = bgsubdirectory+'/'

scaleddirectory = 'scaled'
if not os.path.exists(scaleddirectory):
    print('Creating directory for scaled ota images...')
    os.makedirs(scaleddirectory)

scaledpath = scaleddirectory+'/'

# otastackdirectory = 'otastack'
# if not os.path.exists(otastackdirectory):
#     print 'Creating directory for stacked ota images...'
#     os.makedirs(otastackdirectory)
# 
# otastackpath = otastackdirectory+'/'

skyflatdirectory = 'skyflat'
if not os.path.exists(skyflatdirectory):
    print('Creating directory for sky flats...')
    os.makedirs(skyflatdirectory)

skyflatpath = skyflatdirectory+'/'

coordsdirectory = 'coords'
if not os.path.exists(coordsdirectory):
    print('Creating directory for coordinate files...')
    os.makedirs(coordsdirectory)

coordspath = coordsdirectory+'/'

matchdirectory = 'match'
if not os.path.exists(matchdirectory):
    print('Creating directory for match files...')
    os.makedirs(matchdirectory)

matchpath = matchdirectory+'/'

sdssofflinedir = 'sdssoffline'
if not os.path.exists(sdssofflinedir):
    print('Creating directory for sdss catalogs...')
    os.makedirs(sdssofflinedir)

sdsspath = sdssofflinedir+'/'

# twomassofflinedir = 'twomassoffline'
# if not os.path.exists(twomassofflinedir):
#     print 'Creating directory for 2mass catalogs...'
#     os.makedirs(twomassofflinedir)
# 
# twomasspath = twomassofflinedir+'/'

gaiaofflinedir = 'gaiaoffline'
if not os.path.exists(gaiaofflinedir):
    print('Creating directory for gaia catalogs...')
    os.makedirs(gaiaofflinedir)

gaiapath = gaiaofflinedir+'/'


sourcedir = 'sources'
if not os.path.exists(sourcedir):
    print('Creating directory for detected sources...')
    os.makedirs(sourcedir)

sourcepath = sourcedir+'/'

def cfgparse(cfg_file, verbose=True):
    """
    This function reads the ``yaml`` configuration file that is needed
    to make the ``odi_process.py`` and ``odi_scalestack_process.py`` functions
    work. The configuration file must be given the name ``config.yaml`` and be
    in the same directory as the images that will be processed.

    Parameters
    ----------
    cfg_file : str
        A string containing the name of the configuration file.
    Returns
    -------
    object_str : str
        Name of the object in the images (eg M53)
    filters : list
        List of filters present in images
    instrument : str
        Name of instrument used ``podi`` or ``5odi``
    images : list
        List of images that will be processed
    illcor_flag : bool
        If ``True`` an illumination correction will be applied to the images
    skyflat_src : str
        Can be ``object`` or ``master``. If ``object`` the dark sky flat is
        created from the images being processed. If ``master`` the dark sky
        flat is taken from a master calibration
    wcs_flag : bool
        If ``True`` the a step will be taken to improve the WCS solution on each
        OTA
    reproject_flag : bool
        If ``True`` all of the OTAs will be reprojected to OTA33 of the first
        image in the image list.
    scale_flag : bool
        If ``True`` the OTAs will be scaled to the same level before stacking
    stack_flag : bool
        If ``True`` the processed OTAs will be stacked into a final image
    align_flag : bool
        If ``True`` the stacked images will be aligned using pixel shifts computed from the gaia catalog
    gaia_flag : bool
        If ``True`` the gaia catalog will be used as the sources in the fix WCS
        step
    cluster_flag : bool
        If ``True`` the central regions of a crowded area, such as a globular
        cluster, are avoided when collecting the gaia sources to fix the WCS
    ra_center : float
        If the ``cluster_flag`` is ``True``, then the central ``Ra`` and ``Dec``
        of the cluster, or area that is to be avoided must be given in the
        configuration file.
    dec_center : float
        If the ``cluster_flag`` is ``True``, then the central ``Ra`` and ``Dec``
        of the cluster, or area that is to be avoided must be given in the
        configuration file.
    min_radius : float
        If the ``cluster_flag`` is ``True``, then a minimum radius in
        arcminutes, relative to ``ra_center`` and ``dec_center``, must be given
        in the configuration file.

    """
    from sys import exit
    from yaml import load, dump
    from odi_config import ODIImage
    try:
        from yaml import CLoader as Loader, CDumper as Dumper
    except ImportError:
        from yaml import Loader, Dumper

    with open(cfg_file,'r') as stream:
        data = load(stream, Loader=Loader)
        illcor_flag = (data['processing'])['illumination_correction']
        skyflat_src = (data['processing'])['dark_sky_flat_source']
        scale_flag = (data['processing'])['scale_images']
        wcs_flag = (data['processing'])['wcs_correction']
        reproject_flag = (data['processing'])['reproject']
        stack_flag = (data['processing'])['stack_images']
        align_flag = (data['processing'])['align_images']
        gaia_flag = (data['processing'])['get_gaia']
        cluster_flag = (data['processing'])['cluster_field']
        ra_center = (data['processing'])['ra_center']
        dec_center = (data['processing'])['dec_center']
        min_radius = (data['processing'])['min_radius']


        object_str = (data['basic'])['object']
        filters = (data['basic'])['filters']
        instrument = (data['basic'])['instrument']

        images = {}
        dithers = {}
        scale_ref = {}
        if verbose:
            print('----------------------------------')
            print('odi_tools | Basic data:')
            print('----------------------------------')
            print('object:                 ', object_str)
            print('filters:                ', filters)
            print('instrument:             ', instrument)
            print('----------------------------------')
            print('Steps to perform:')
            print('----------------------------------')
            print('illumination correction:', illcor_flag)
            print('dark sky flat source:   ', skyflat_src)
            print('wcs correction:         ', wcs_flag)
            print('reprojection:           ', reproject_flag)
            print('scaling:                ', scale_flag)
            print('stacking:               ', stack_flag)
            print('aligning:               ', align_flag)
            print('----------------------------------')
            print('Images (* = scaling ref. image)')
            print('----------------------------------')
        header_string = 'dither  '
        for filter in filters:
            imglist = {}
            try:
                for d,f in data[filter].items():
                    # only keep integer numbers as dither names, 
                    if type(d) is int:
                        im = ODIImage(f, d, instrument)
                        imglist[d] = im
                    # keep the 'ref' in a dictionary for passing. other string names are ignored!
                    elif 'ref' in d:
                        scale_ref[filter] = imglist[f]
                images[filter] = list(imglist.values())
            except KeyError:
                print("images for filter '"+filter+"' not defined in configuration file...")
                exit()
            header_string = header_string + filter + ' '+' '*(len(data[filter][1])-len(filter)+1)
            # print scale_ref[filter]
        if verbose:
            print(header_string)
            dithernos = set()
            for filt in filters:
                dithernos = dithernos | set(data[filt].keys())
            # remove the 'ref' designation from the set of dither numbers if it exists
            if 'ref' in dithernos:
                dithernos.remove('ref')
            for dither in dithernos:
                dither_string = '   {:2d}  '.format(dither)
                for filter in filters:
                    try:
                        if filter in list(scale_ref.keys()) and data[filter][dither] == scale_ref[filter].f:
                            dither_string = dither_string + '*'+data[filter][dither]+' '
                        else:
                            dither_string = dither_string + ' '+data[filter][dither]+' '
                    except KeyError:
                        dither_string = dither_string + ' --no data'+'-'*(len(data[filter][1])-9)+' '
                print(dither_string)        
        return object_str, filters, instrument, images, illcor_flag, skyflat_src, wcs_flag, reproject_flag, scale_flag, scale_ref, stack_flag, align_flag, gaia_flag, cluster_flag, ra_center, dec_center, min_radius

def photcfgparse(cfg_file):
    """
    This function reads the ``yaml`` configuration file that is needed
    to make the ``odi_phot_process.py`` function work. The configuration file
    must be given the name ``phot_config.yaml`` and be in the same directory
    as the images that will be processed.

    Parameters
    ----------
    cfg_file : str
        A string containing the name of the configuration file.
    Returns
    -------
    object_str : str
        Name of the object in the images (eg M53)
    filters : list
        List of filters present in images
    instrument : str
        Name of instrument used ``podi`` or ``5odi``
    images : list
        List of images that will be processed
    new_extension : str
        The new extension that will be given to the images resulting from
        ``odi_phot_process.py``
    remove_tpv_flag : bool
        If ``True`` header cards with ``TPV`` will be removed from the stacked
        images previously produced by ``odi_process.py``
    trim_image_flag : bool
        If ``True`` the stacked images previously produced by
        ``odi_process.py`` will be trimmed.
    wcs_flag : bool
        If ``True`` the a step will be taken to improve the WCS on
        the stacked images previously produced by
        ``odi_process.py``
    trim_image_flag : bool
        The sections to trim from the stacked images.
    airmasses : list
        The airmasses of the images in the order they will be processed.
    """
    from sys import exit
    from yaml import load, dump
    try:
        from yaml import CLoader as Loader, CDumper as Dumper
    except ImportError:
        from yaml import Loader, Dumper
    with open(cfg_file,'r') as stream:
        data = load(stream, Loader=Loader)
        remove_tpv_flag = (data['processing'])['remove_tpv']
        trim_image_flag = (data['processing'])['trim_image']
        wcs_flag = (data['processing'])['wcs_correction']
        trim_section = (data['processing'])['trim_section']
        airmasses = (data['processing'])['airmasses']
        new_extension = (data['processing'])['new_extension']

        object_str = (data['basic'])['object']
        filters = (data['basic'])['filters']
        instrument = (data['basic'])['instrument']

        print('----------------------------------')
        print('odi_tools | Basic data:')
        print('----------------------------------')
        print('object:                 ', object_str)
        print('filters:                ', filters)
        print('instrument:             ', instrument)
        print('trim section:           ', trim_section)
        print('arimasses:              ', airmasses)
        print('----------------------------------')
        print('Steps to perform:')
        print('----------------------------------')
        print('remove TPV header keys:', remove_tpv_flag)
        print('fix wcs solution:      ', wcs_flag)
        print('trim image:            ', trim_image_flag)
        print('new extension:         ', new_extension)
        print('----------------------------------')
        print('Images:')
        print('----------------------------------')

        images = {}
        for filter in filters:
            try:
                images[filter] = data[filter]
                # print data[filter].keys()
            except KeyError:
                print("images for filter '"+filter+"' not defined in configuration file...")
                exit()
            print(images[filter][1])

    return object_str, filters, instrument, images, new_extension, remove_tpv_flag, trim_image_flag, wcs_flag, trim_section, airmasses

def main():
    object_str, filters, instrument, images, illcor_flag, skyflat_src, wcs_flag, reproject_flag, scale_flag, scale_ref, stack_flag, align_flag, gaia_flag, cluster_flag, ra_center, dec_center, min_radius = cfgparse('example_config.yaml', verbose=True)
    print(scale_ref)

if __name__ == '__main__':
    main()
