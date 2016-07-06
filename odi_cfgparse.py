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
        reproj_flag = (data['processing'])['reproject']
        scale_flag = (data['processing'])['scale_images']
        stack_flag = (data['processing'])['stack_images']
        
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
        print 'reprojection:           ', reproj_flag
        print 'scaling:                ', scale_flag
        print 'stacking:               ', stack_flag
        print '----------------------------------'
        print 'Images:'
        print '----------------------------------'
        header_string = 'dither '
        for filter in filters:
            header_string = header_string + filter + '                                   '
            try:
                images[filter] = data[filter]
                # print data[filter].keys()
            except KeyError:
                print "images for filter '"+filter+"' not defined in configuration file..."
                exit()
        print header_string
        for dither in images[filters[0]].keys():
            dither_string = '     '+repr(dither)+' '
            for filter in filters:
                dither_string = dither_string + images[filter][dither]+' '
            print dither_string
        return object_str, filters, instrument, images, illcor_flag, skyflat_src, reproj_flag, scale_flag, stack_flag
        

def main():
    pass

if __name__ == '__main__':
    cfgparse('default_config.yaml')