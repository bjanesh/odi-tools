.. _basic_usage:

Basic usage
===========

All you need to do to get started is download your QR-ed data from the ODI-PPA
using the wget download command, then follow these steps. An explination of
running quick reduce from ODI-PPA will be given in other sections of the
documentation:

*******************
Preparing your data
*******************
1. move all individual ``.fz`` files into the top level folder:
   ``mv calibrated/**/*.fz .``
2. unpack the compressed fits files using `funpack
   <https://heasarc.gsfc.nasa.gov/fitsio/fpack/>`_
3. you need to rename your files to match the appropriate dither
   pointing identification. for example, QR files are named by the pattern
   ``OBSID_OBJECT_FILTER.JOBID.fits``. The final digit of the OBSID e.g.
   ``20151008T195949.1`` needs to match the number sequence of the dithers 1-9.
   Your data may not match this pattern due to restarted observations, multiple
   night observations, etc.

***********************************
Running the code (a broad overview)
***********************************

1. copy ``example_config.yaml`` to your data directory as ``config.yaml`` and edit the
   file to match your preferences/data. Make sure that the number for each image
   matches the correct number in the dither sequence!
2. run ``odi_process.py`` in the folder containing the unpacked/renamed fits images.
   This will (optionally) illumination correct the images, fix their WCS,
   reproject them to a common pixel scale, and perform background subtraction on them.
3. this will take a while, so make sure nothing bad happened
4. run ``odi_scalestack_process.py`` in the folder containing the unpacked/renamed
   fits images. This will detect bright stellar sources in the images and use
   them to calculate a scaling factor relative to the image in the sequence
   with the lowest airmass, then apply the scale, stack the images,
   then add in a common background value.
5. finished! check your images to make sure everything went okay.

.. _example_config:

**************************
Example configuration file
**************************

Here are the contents of ``example_config.yaml`` available on the `odi-tools
GitHub repo
<https://github.com/bjanesh/odi-tools/blob/master/example_config.yaml>`_
::

    # odi-tools configuration file
    basic:
      object: M13                           # the name of your object
      filters: [odi_g, odi_r, odi_i]        # correct filter strings required
      instrument: 5odi                      # podi | 5odi | mosaic; script will
                                            # verify using image header info

    processing:                             # optional steps performed in odi_process.py
      illumination_correction: yes          # if yes, set dark sky flat source below
      dark_sky_flat_source: object          # object | master
      wcs_correction: yes
      reproject: yes
      scale_images: yes
      stack_images: yes

    # list the images you want to process
    # be sure to associate the filename with the correct dither pointing!
    # OBSID and image header are NOT always an accurate reflection of the absolute dither position
    # so you must use your notes / observing log to define them here
    # sections must be named according to the filter names

    odi_g:
      1: 20130510T002928.1_m13-9_odi_g.5869.fits
      2: 20130510T002928.2_m13-9_odi_g.5869.fits
      3: 20130510T002928.3_m13-9_odi_g.5869.fits
      4: 20130510T002928.4_m13-9_odi_g.5869.fits
      5: 20130510T002928.5_m13-9_odi_g.5869.fits
      6: 20130510T002928.6_m13-9_odi_g.5869.fits
      7: 20130510T002928.7_m13-9_odi_g.5869.fits
      8: 20130510T002928.8_m13-9_odi_g.5869.fits
      9: 20130510T002928.9_m13-9_odi_g.5869.fits

    odi_r:
      1: 20130510T002928.1_m13-9_odi_r.5869.fits
      2: 20130510T002928.2_m13-9_odi_r.5869.fits
      3: 20130510T002928.3_m13-9_odi_r.5869.fits
      4: 20130510T002928.4_m13-9_odi_r.5869.fits
      5: 20130510T002928.5_m13-9_odi_r.5869.fits
      6: 20130510T002928.6_m13-9_odi_r.5869.fits
      7: 20130510T002928.7_m13-9_odi_r.5869.fits
      8: 20130510T002928.8_m13-9_odi_r.5869.fits
      9: 20130510T002928.9_m13-9_odi_r.5869.fits

    odi_i:
      1: 20130510T002928.1_m13-9_odi_i.5869.fits
      2: 20130510T002928.2_m13-9_odi_i.5869.fits
      3: 20130510T002928.3_m13-9_odi_i.5869.fits
      4: 20130510T002928.4_m13-9_odi_i.5869.fits
      5: 20130510T002928.5_m13-9_odi_i.5869.fits
      6: 20130510T002928.6_m13-9_odi_i.5869.fits
      7: 20130510T002928.7_m13-9_odi_i.5869.fits
      8: 20130510T002928.8_m13-9_odi_i.5869.fits
      9: 20130510T002928.9_m13-9_odi_i.5869.fits
