from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.visualization import *
from astropy.visualization.mpl_normalize import ImageNormalize

def wcs_header(image,ota):
  """
  A function to update the keywords in a QR'd image that start
  with PV to TPV. This is done in order to be used with the
  astropy wcs and reproject packages.
  
  image: fits image to be used
  ota: ota you would like to update. The qr'd images come down
  as a multi extension fits file. extenstions 1-9 are the 9
  different ota's. 
  
  This fucntion will return a dictionary 'header'. This will
  not overwrite the header for any given ota, but can be passed
  to the various atropy packages
  
  """
  hdu = fits.open(image)
  header = hdu[ota].header
  oldkeys = []
  newkeys = []

  for key,i in header.iteritems():
    if 'PV' in key[0:2]:
      oldkeys.append(key)
      newkeys.append('T'+key)

  for i in range(len(oldkeys)):
    header[newkeys[i]] = header.pop(oldkeys[i])
    
  return hdu,header

# Example use, imshow will have ra dec in degrees rather than x and y.    
hdu,header = wcs_header('20140406T225633.1_GCPair-F1_odi_r.5898.fits.fz',1)
norm = ImageNormalize(stretch=LinearStretch())
ax1 = plt.subplot(1,1,1, projection=WCS(header))
ax1.imshow(hdu[1].data,cmap='Greys_r',vmin=850,vmax=1065,origin='lower',norm=norm)
plt.show()


  
