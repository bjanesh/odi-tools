#!/usr/bin/env python
import os
import sys
import shutil
import glob
import odi_config as odi

def main():
    try:
        object_str, filters, instrument, images, illcor_flag, skyflat_src, wcs_flag, reproject_flag, scale_flag, stack_flag, gaia_flag, cluster_flag, ra_center, dec_center, min_radius = odi.cfgparse('config.yaml')
    except IOError:
        print 'config.yaml does not exist, quitting...'
        exit()
    object_str = object_str.replace(' ','')
    # os.mkdir(object_str)
#     shutil.copy('derived_props.txt', object_str)
#shutil.copy('config.yaml', object_str)
    fits_files = glob.glob(object_str+'*.fits')
    pl_files = glob.glob(object_str+'*_bpm.pl')
    scale_files = glob.glob('*scales.txt')
    print fits_files, pl_files, scale_files
    pass

if __name__ == '__main__':
    main()
