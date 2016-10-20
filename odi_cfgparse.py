def cfgparse(cfg_file):
    from sys import exit
    from yaml import load, dump
    try:
        from yaml import CLoader as Loader, CDumper as Dumper
    except ImportError:
        from yaml import Loader, Dumper

    with open(cfg_file,'r') as stream:
        data = load(stream, Loader=Loader)
        illcor_flag = (data['processing'])['illumination_correction']
        skyflat_src = (data['processing'])['dark_sky_flat_source']
        scale_flag = (data['processing'])['scale_images']
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
        print 'scaling:                ', scale_flag
        print 'stacking:               ', stack_flag
        print '----------------------------------'
        print 'Images:'
        print '----------------------------------'
        header_string = 'dither '
        for filter in filters:
            try:
                images[filter] = data[filter]
                # print data[filter].keys()
            except KeyError:
                print "images for filter '"+filter+"' not defined in configuration file..."
                exit()
            header_string = header_string + filter + ' '*(len(images[filter][1])-len(filter)+1)
        print header_string
        for dither in images[filters[0]].keys():
            dither_string = '     '+repr(dither)+' '
            for filter in filters:
                dither_string = dither_string + images[filter][dither]+' '
            print dither_string
        return object_str, filters, instrument, images, illcor_flag, skyflat_src, scale_flag, stack_flag, gaia_flag, cluster_flag, ra_center, dec_center, min_radius

def photcfgparse(cfg_file):
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
    cfgparse('default_config.yaml')
