.. _odi_process_bd:

=================================
A detailed look at odi_process.py
=================================

``odi_process.py`` is the script responsible for carrying out all of
the steps in the ``odi-tools`` pipeline. Here we will give a detailed
explanation about each step the the script is doing.

Reading the configuration file
------------------------------
The first thing ``odi_process.py`` does is try to
read and parse the configuration file that should
also be located in the current working directory. This
is done with the function ``odi.cfgparse``.
This file has to be called ``config.yaml``.
If this file is not found, the program will exit and
the user should ensure their configuration file is present
and been given the right name. These line are responsible
for creating variables that will be needed for the
rest of the pipeline to function.

Here is a list of the variables set by the configuration file.

- object_str
- filters
- instrument
- images
- illcor_flag
- skyflat_src
- wcs_flag
- reproject_flag
- scale_flag
- stack_flag
- gaia_flag
- cluster_flag
- ra_center
- dec_center
- min_radius

Creating the image lists
------------------------
The next step in ``odi_process.py`` is to create
the list of images that will be processed. This list
is given then name ``images_``. The list is populated
by iterating the ``images`` dictionary returned
by the previous ``odi.cfgparse`` step. The
items in ``images_`` will be in the same order as
they appear in ``config.yaml`` and separated by filter.

Setting the reprojection image, source catalog, and instrument
--------------------------------------------------------------
All of the images processed ``odi_process.py`` will
be reprojected according to OTA33 in the first image
in the ``images_`` list. This should correspond to the
first image in your dither pattern for the set of
images your are currently processing. The coordinates of
this OTA in this image are ``rad`` and ``decd`` and they
are returned by the function ``odi.get_targ_ra_dec``.

In order to improve the WCS solution on each OTA, ``odi_process.py``
requires a source catalog with known Ra and Dec values. To set the
desired source catalog, ``odi_process.py`` checks if the user has
has set the ``gaia_flag`` to ``yes`` in ``config.yaml``. If this is
the case then ``odi_process.py`` will use the Gaia catalog as
the source list for the fixing the WCS solutions. If the
``gaia_flag`` is set to ``no``, ``odi_process.py`` will default
to using the SDSS catalog. This step sets the ``source`` variable.

For ``odi_process.py`` to run correctly, the pipeline must
also be told if the data being processed are from ``pODI`` or
``ODI``. This is accomplished by the ``odi.instrument`` function
that reads the header of the first item in the ``images_`` list
and returns the ``inst`` variable.

Creating source catalogs
------------------------
