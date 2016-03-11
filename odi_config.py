from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.io import fits
from astropy.wcs import WCS
from rand_bkg import bkg_boxes
from astropy.convolution import Gaussian2DKernel
from astropy.stats import sigma_clipped_stats
from photutils.detection import detect_sources
from photutils.utils import random_cmap
from scipy.ndimage import binary_dilation

from odi_calibrate import *
from odi_helpers import *
from mask_ota import *
from get_gaps import *
from bpm_tools import *
from fixwcs import *
from illcor import *
from getfwhm import *
from zeropoint import *
from scale import *
from sdss_coords import *
from new_scale import *
from offlinecats import *
from ota_sourcefind import *

OTA_dictionary = {1:'OTA33.SCI',2: 'OTA34.SCI',3 :'OTA44.SCI', 4:'OTA43.SCI',5:'OTA42.SCI', 6:'OTA32.SCI', 
7:'OTA22.SCI' ,8:'OTA23.SCI',9:'OTA24.SCI'}

bpmdirectory = 'bpmasks'
if not os.path.exists(bpmdirectory):
	print 'Creating directory for bad pixel masks...'
	os.makedirs(bpmdirectory)

bppath = bpmdirectory+'/'

illcordirectory = 'illcor'
if not os.path.exists(illcordirectory):
	print 'Creating directory for illumination corrected ota images...'
	os.makedirs(illcordirectory)

illcorpath = illcordirectory+'/'

reprojdirectory = 'reproj'
if not os.path.exists(reprojdirectory):
	print 'Creating directory for reprojected ota images...'
	os.makedirs(reprojdirectory)

reprojpath = reprojdirectory+'/'

bgsubdirectory = 'bgsub'
if not os.path.exists(bgsubdirectory):
	print 'Creating directory for background subtracted ota images...'
	os.makedirs(bgsubdirectory)

bgsubpath = bgsubdirectory+'/'

scaleddirectory = 'scaled'
if not os.path.exists(scaleddirectory):
	print 'Creating directory for scaled ota images...'
	os.makedirs(scaleddirectory)

scaledpath = scaleddirectory+'/'

otastackdirectory = 'otastack'
if not os.path.exists(otastackdirectory):
	print 'Creating directory for stacked ota images...'
	os.makedirs(otastackdirectory)

otastackpath = otastackdirectory+'/'

skyflatdirectory = 'skyflat'
if not os.path.exists(skyflatdirectory):
	print 'Creating directory for sky flats...'
	os.makedirs(skyflatdirectory)

skyflatpath = skyflatdirectory+'/'

coordsdirectory = 'coords'
if not os.path.exists(coordsdirectory):
	print 'Creating directory for coordinate files...'
	os.makedirs(coordsdirectory)

coordspath = coordsdirectory+'/'

matchdirectory = 'match'
if not os.path.exists(matchdirectory):
	print 'Creating directory for match files...'
	os.makedirs(matchdirectory)

matchpath = matchdirectory+'/'

sdssofflinedir = 'sdssoffline'
if not os.path.exists(sdssofflinedir):
	print 'Creating directory for sdss catalogs...'
	os.makedirs(sdssofflinedir)

sdsspath = sdssofflinedir+'/'

twomassofflinedir = 'twomassoffline'
if not os.path.exists(twomassofflinedir):
	print 'Creating directory for 2mass catalogs...'
	os.makedirs(twomassofflinedir)

twomasspath = twomassofflinedir+'/'


sourcedir = 'sources'
if not os.path.exists(sourcedir):
	print 'Creating directory for detected sources...'
	os.makedirs(sourcedir)

sourcepath = sourcedir+'/'