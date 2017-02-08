#!/usr/bin/env python
import os
import sys
import pwd
import shutil
import glob
import datetime
import tarfile
import odi_config as odi
from subprocess import call

def main():
    try:
        object_str, filters, instrument, images, illcor_flag, skyflat_src, wcs_flag, reproject_flag, scale_flag, stack_flag, gaia_flag, cluster_flag, ra_center, dec_center, min_radius = odi.cfgparse('config.yaml', verbose=False)
    except IOError:
        print 'config.yaml does not exist, quitting...'
        exit()
    
    funpack_path = '/usr/bin/fpack'
    
    # replace any spaces in the object name with -nothing-
    object_str = object_str.replace(' ','')
    
    # get today's date
    today = datetime.date.today()
    
    # get the username of the owner of the folder
    stat_info = os.stat('.')
    uid = stat_info.st_uid
    user = pwd.getpwuid(uid)[0]
    
    # construct the filename --> yyyy-mm-dd_username_object
    file_stem =  '{:0>4}-{:0>2}-{:0>2}'.format(today.year,today.month,today.day)+'_'+user+'_'+object_str
    
    # figure out the fits, pl, scale filenames (1/filter)
    fits_files = glob.glob(object_str+'*.fits')
    pl_files = glob.glob(object_str+'*_bpm.pl')
    scale_files = glob.glob('*scales.txt')
    
    # use fpack to compress the fits files in a sane way
    for file_ in fits_files:
        if not os.path.isfile(file_+'.fz'):
            funpack_cmd = funpack_path+' '+file_+'.fz'
            call(funpack_cmd, shell=True)
    
    # if there's not already an archive here...
    if not os.path.isfile('/ssd1/'+file_stem+".tar"):
        # make a tar object to move the kept files into
        tar = tarfile.open('/ssd1/'+file_stem+".tar","w")
    
        # move the derived props and config files
        tar.add('derived_props.txt')
        tar.add('config.yaml')
    
        # move the fits files
        for f in fits_files:
            tar.add(f+'.fz')
    
        # move the pl files    
        for p in pl_files:
            tar.add(p)
    
        # move the scale files
        for s in scale_files:
            tar.add(s)
        
        tar.close()
        
    print 'Files contained in /ssd1/'+file_stem+".tar:"
    print '- config.yaml'
    print '- derived_props.txt'
    # move the fits files
    for f in fits_files:
        print '- '+f

    # move the pl files    
    for p in pl_files:
        print '- '+p

    # move the scale files
    for s in scale_files:
        print '- '+s
    print 'To copy this to your computer, run:'
    print '--> scp /ssd1/'+file_stem+".tar USER@HOST.astro.indiana.edu:~"
    print '  where USER and HOST are your local user name and computer name.'
    print '  This will place the file in your local home directory.'
    print '  You can change this directory if you like (e.g. ~/images).'
    print 'Then notify Bob Lezotte (hlezotte [at] indiana.edu)'
    print '  that you are finished processing with wopr.'
    
    pass

if __name__ == '__main__':
    main()
