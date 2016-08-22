Navigating a QR image
=====================


.. _opening-the-fits-file:

Opening the fits file
---------------------

``odi-tools`` utilizes the ``fits.io`` module in the ``astropy`` package to
open the multi-extension QR fits files:

>>> from astropy.io import fits
>>> img = '20140406T214040.2_GCPair-F1_odi_g.6183.fits'
>>> hdulist = fits.open(img)


.. _accessing-the-extensions:

Accessing the extensions
------------------------

Looking at ``hdulist`` is now an ``astropy`` HDUlist where each element
cooresponds to an extension of the fits file.

>>> print hdulist.info()
Filename: 20140406T214040.2_GCPair-F1_odi_g.6183.fits
No.    Name         Type      Cards   Dimensions     Format
0    PRIMARY     PrimaryHDU     294   ()
1    OTA33.SCI   ImageHDU       285   (4096, 4096)   float32
2    OTA34.SCI   ImageHDU       285   (4096, 4096)   float32
3    OTA44.SCI   ImageHDU       285   (4096, 4096)   float32
4    OTA43.SCI   ImageHDU       285   (4096, 4096)   float32
5    OTA42.SCI   ImageHDU       285   (4096, 4096)   float32
6    OTA32.SCI   ImageHDU       285   (4096, 4096)   float32
7    OTA22.SCI   ImageHDU       285   (4096, 4096)   float32
8    OTA23.SCI   ImageHDU       285   (4096, 4096)   float32
9    OTA24.SCI   ImageHDU       321   (4096, 4096)   float32
10   CAT.2MASS   BinTableHDU     15   680R x 2C    [D, D]
11   CAT.ODI     BinTableHDU    138   562R x 33C   [D, D, E, E, E, E, E, I, I, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E]
12   CAT.ODI+2MASS  BinTableHDU  31   165R x 5C    [D, D, D, D, E]
13   CAT.PHOTCALIB  BinTableHDU 172   251R x 36C   [D, D, D, D, E, E, E, E, E, E, E, E, E, E, D, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E]
14   SKYLEVEL    BinTableHDU     50   1663R x 7C   [D, D, D, D, D, D, I]
15   ASSOCIATIONS  BinTableHDU   17   16R x 3C     [25A, 375A, 100A]

This particular image  is from ``podi`` meaning there are a total of 9 OTA in
the focal plane. In the ``hdulist`` the extensions for each of the OTAs are
given in rows 1-9 and have names with the following convention ``OTAxy.SCI``.
For a 5x6 ``ODI`` image, there would be 30 of these extensions.

The ``hdulist`` allows us to access each of the extensions by the its ``name``
give in the Name column, or by its number. If, for example, we wanted to only
look at ``OTA33`` we can do the following:

>>> hdu_ota = hdulist['OTA33.SCI']
# or using the number
>>> hdu_ota = hdulist[1]

Now that we have isolated this single OTA we can pick out its individual
header and data:

>>> ota_header = hdu_ota.header
>>> ota_data = hdu_ota.data

``ota_header`` and ``ota_data`` are now easily accessible to use in other
scripts and modules. It is also easy to access ``header`` keywords:

>>> print ota_header['CRVAL1']
198.750082573
>>> print ota_header['CRVAL2']
18.4113084341

It is also convenient to pass individual OTAs to ``pyraf`` tasks based on the
name of the extension. If we wanted to run ``daofind`` on OTA33 we could do
the following:

>>> from pyraf import iraf
>>> img = '20140406T214040.2_GCPair-F1_odi_g.6183.fits'
>>> iraf.unlearn(iraf.apphot.daofind)
>>> iraf.datapars.setParam('fwhmpsf',fwhm,check=1)
>>> iraf.datapars.setParam('datamin',-900,check=1)
>>> iraf.datapars.setParam('datamax',60000,check=1)
>>> iraf.datapars.setParam('sigma',bg_std,check=1)
>>> iraf.findpars.setParam('threshold',2.5)
>>> iraf.apphot.daofind.setParam('output',output.txt)
>>> iraf.apphot.daofind(image=img+'['+'OTA33.SCI'+']', verbose="no", verify='no')

.. _odi-tools-scheme:

The odi-tools scheme
--------------------

The modules in ``odi-tools`` are designed to work over lists of images while
stepping through each of the OTA extensions in a given image. This will
discussed in further detail in other parts of the documentation. An pseudo
code of this scheme would be:

>>> imglist = ['img1.fits','img2.fits','img3.fits']
>>> ota_dictionary = {1:'OTA33.SCI',2: 'OTA34.SCI',3 :'OTA44.SCI',
...                   4:'OTA43.SCI',5:'OTA42.SCI', 6:'OTA32.SCI',
...                   7:'OTA22.SCI' ,8:'OTA23.SCI',9:'OTA24.SCI'}
>>> for img in imglist:
...     for key in ota_dictionary:
...         ota = ota_dictionary[key]
...         perform tasks on img[ota]

.. _the-other-extensions:

The other extensions
--------------------

In addition to the extensions for each OTA, the ``hdulist`` also contains
extensions linking to fits tables with useful information. They are
``CAT.2MASS``, ``CAT.ODI``, ``CAT.ODI+2MASS``, ``CAT.PHOTCATLIB``, ``SKYLEVEL``,
``ASSOCIATIONS``. The header and data in each of these tables are easily accessed.

>>> photcat_data = hdulist['CAT.PHOTCATLIB'].data
>>> photcat_header = hdulist['CAT.PHOTCATLIB'].header

Some of the information in these tables is used y ``odi-tools`` during the
image processing.
