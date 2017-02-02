Python configuration to run odi-tools
=====================================

The function ``odi_config.py`` imports all of the modules needed to run
``odi-tools`` as well as sets up the needed dictionaries and directories.
It is imported in the following way in the main ``odi-tool`` scripts
(``odi_process``, ``odi_scalestack_process``, ``odi_phot_process``).

>>> import odi_config as odi

Once it is imported, the dictionaries are directories created in
``odi_config.py`` can be referenced throughout the ``odi-tools`` pipeline in
the following manner

>>> example_dict = odi.dictionary
>>> example_directory = odi.directory

Similarly, we can also reference all of the modules and functions imported
in ``odi_config.py``.

>>> gaps = odi.get_gaps(img, ota)

The OTA dictionaries
--------------------
There are two dictionaries defined by this function that correspond to
different versions of ODI.

1. ``podi_dictionary``
2. ``odi5_dictionary``

As an example, here are the contents of ``podi_dictionary``::

    podi_dictionary = {1: 'OTA33.SCI',
                       2: 'OTA34.SCI',
                       3: 'OTA44.SCI',
                       4: 'OTA43.SCI',
                       5: 'OTA42.SCI',
                       6: 'OTA32.SCI',
                       7: 'OTA22.SCI',
                       8: 'OTA23.SCI',
                       9: 'OTA24.SCI'
                       }

For those that are not familiar with Python, a dictionary is made up of pairs
of ``keys`` and ``values``. The ``keys`` are to the left of the colon, and the
``values`` are to the right. In our case, the ``keys`` are the numbers 1-9, and
the values are the names of the different OTAs, (e.g. ``OTA33.SCI``). The
dictionaries provide a clean way to work through a multi-extension fits image
like those produced by ODI. The simpled code example below provides an
illustration of how this is done in ``odi-tools``.

>>> images = ['img1.fits','img2.fits']
>>> for img in images:
>>>     for key in podi_dictionary:
>>>         print img, key

This would produce the following output::

'img1.fits' 'OTA33.SCI'
'img1.fits' 'OTA34.SCI'
'img1.fits' 'OTA44.SCI'
'img1.fits' 'OTA43.SCI'
'img1.fits' 'OTA42.SCI'
'img1.fits' 'OTA32.SCI'
'img1.fits' 'OTA22.SCI'
'img1.fits' 'OTA23.SCI'
'img1.fits' 'OTA24.SCI'
'img2.fits' 'OTA33.SCI'
'img2.fits' 'OTA34.SCI'
'img2.fits' 'OTA44.SCI'
'img2.fits' 'OTA43.SCI'
'img2.fits' 'OTA42.SCI'
'img2.fits' 'OTA32.SCI'
'img2.fits' 'OTA22.SCI'
'img2.fits' 'OTA23.SCI'
'img2.fits' 'OTA24.SCI'

Although this is a simple example it illustrates the overall workflow of
``odi-tools``.

The ``odi5_dictionary`` works the same way, but simply has more OTAs. The
correct dictionary is selected by :py:func:`odi_helpers.instrument`. If the
instrument used was ``5odi``, ``odi5_dictionary`` is used for ``odi-tools``,
if it was ``podi``, ``podi_dictionary`` is used.

Processing directories
----------------------

``odi_config.py`` also sets up a number of directories to hold the
intermediate data products during the data processing. Here is an example
of how one of those directories is created

>>> bpmdirectory = 'bpmasks'
>>> if not os.path.exists(bpmdirectory):
>>>    print 'Creating directory for bad pixel masks...'
>>>    os.makedirs(bpmdirectory)

>>> bppath = bpmdirectory+'/'

The directory in this case is given the name ``bpmasks``. Then, we check if
the directory already exists. If is does not, the directory is created. Once
this directory is created it can be accessed by other ``odi-tools`` modules
and scripts using the following

>>> odi.bppath

Here is a full list of the directories created

- ``bpmasks`` - directory for bad pixel masks
- ``illcor`` - directory for illumination corrected ota images
- ``reproj`` - directory for reprojected ota images
- ``bgsub`` - directory for background subtracted ota images
- ``scaled`` - directory for scaled ota images
-  ``otastack`` - directory for stacked ota images
- ``skyflat`` - directory for sky flats
- ``coords`` - directory for coordinate files
- ``match`` - directory for match files
- ``sdssoffline``- directory for sdss catalogs
- ``twomassoffline`` - directory for 2mass catalogs
- ``gaiaoffline`` - directory for gaia catalogs
- ``sources`` - directory for detected sources
