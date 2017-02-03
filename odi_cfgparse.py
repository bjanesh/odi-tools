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
        if verbose:
            print '----------------------------------'
            print 'odi_tools | Basic data:'
            print '----------------------------------'
            print 'object:                 ', object_str
            print 'filters:                ', filters
            print 'instrument:             ', instrument
            print '----------------------------------'
            print 'Steps to perform:'
            print '----------------------------------'
            print 'illumination correction:', illcor_flag
            print 'dark sky flat source:   ', skyflat_src
            print 'wcs correction:         ', wcs_flag
            print 'reprojection:           ', reproject_flag
            print 'scaling:                ', scale_flag
            print 'stacking:               ', stack_flag
            print '----------------------------------'
            print 'Images:'
            print '----------------------------------'
        header_string = 'dither '
        for filter in filters:
            imglist = []
            try:
                for d,f in data[filter].iteritems():
                    imglist.append(ODIImage(f, d, instrument))
                images[filter] = imglist
            except KeyError:
                print "images for filter '"+filter+"' not defined in configuration file..."
                exit()
            header_string = header_string + filter + ' '*(len(data[filter][1])-len(filter)+1)
        if verbose:
            print header_string
            dithernos = set()
            for filt in filters:
                dithernos = dithernos | set(data[filt].keys())
            for dither in dithernos:
                dither_string = '     '+repr(dither)+' '
                for filter in filters:
                    try:
                        dither_string = dither_string + data[filter][dither]+' '
                    except KeyError:
                        dither_string = dither_string + '--no data'+'-'*(len(data[filter][1])-9)+' '
                print dither_string        
        return object_str, filters, instrument, images, illcor_flag, skyflat_src, wcs_flag, reproject_flag, scale_flag, stack_flag, gaia_flag, cluster_flag, ra_center, dec_center, min_radius

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

        print '----------------------------------'
        print 'odi_tools | Basic data:'
        print '----------------------------------'
        print 'object:                 ', object_str
        print 'filters:                ', filters
        print 'instrument:             ', instrument
        print 'trim section:           ', trim_section
        print 'arimasses:              ', airmasses
        print '----------------------------------'
        print 'Steps to perform:'
        print '----------------------------------'
        print 'remove TPV header keys:', remove_tpv_flag
        print 'fix wcs solution:      ', wcs_flag
        print 'trim image:            ', trim_image_flag
        print 'new extension:         ', new_extension
        print '----------------------------------'
        print 'Images:'
        print '----------------------------------'

        images = {}
        for filter in filters:
            try:
                images[filter] = data[filter]
                # print data[filter].keys()
            except KeyError:
                print "images for filter '"+filter+"' not defined in configuration file..."
                exit()
            print images[filter][1]

    return object_str, filters, instrument, images, new_extension, remove_tpv_flag, trim_image_flag, wcs_flag, trim_section, airmasses

def main():
    pass

if __name__ == '__main__':
    cfgparse('example_config.yaml')
