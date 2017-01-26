---
---
### Processing details
*last modified 2 May 2016*

1. Add the images you wish to reduce to a collection in the ODI-PPA and from the collection action menu, choose QuickReduce.
2. Run QuickReduce with the following options: <br />☑️ WCS <br />☑️ Photometry <br />☑️ Fringe (i- and z-band only) <br />☑️ Persistency <br />☑️ Nonlinearity <br />☑️ Pupil Ghost (calibrations) <br />❌ Pupil Ghost (science) <br />☑️ Use Bad Pixel Masks <br />☑️ Cosmic Ray Removal, 3 iteration(s) <br />Make sure you name the job and check the "Email me when finished" box. Wait.
3. Download the QuickReduce data to your hard drive using the "Download Results" button on the QuickReduce job page. You don't need to select any of the optional boxes, just name the job and click submit. Eventually a `wget` command will pop up. Copy it to your clipboard, navigate to the folder you want the data to go into, then paste the `wget` command in your command line.
4. Move all individual .fz files into the top level folder: `mv calibrated/**/*.fz .`
5. Unpack the compressed fits files using `funpack` [link](https://heasarc.gsfc.nasa.gov/fitsio/fpack/)
6. You need to rename your files to match the appropriate dither pointing identification. for example, QR files are named by the pattern `OBSID_OBJECT_FILTER.JOBID.fits`. The final digit of the OBSID `e.g. 20151008T195949.1` needs to match the number sequence of the dithers 1-9. _Your data may not match this pattern_ due to restarted observations, multiple night observations, etc. Support will be added in the future for FITS header keywords providing this functionality without user interaction.
7. Once your files are uncompressed and renamed, run `odi_process.py` in the folder with the images. This script will do the following (most of these steps are applied on an OTA by OTA basis, not the entire image at once):
  1. create catalogs of SDSS and 2MASS stars from the QuickReduce header tables. These will be used later for a bunch of things.
  2. create lists of all object images
  3. open each object image and create a mask of bad pixels _and a mask of all objects_ in the image, write the name of these masks to the header.
  4. take the object images and the bad pixel/object masks and do a masked median combine to create a dark sky flat (one per filter).
  5. illumination correction -- divide each object image by its filter's dark sky flat.
  6. use the SDSS coordinates from before to select bright stars and use them to correct the WCS in each image using `iraf.mscred.msccmatch` (this works poorly in crowded fields at the moment). Do three iterations of `fix_wcs`.
  7. now that the WCS is good in each image, reproject each image of the object in all filters to a common system, using the _target_ RA and Dec as the reference point. (using `iraf.mscred.mscimage`)
  8. measure the mean stellar fwhm in each OTA using `iraf.imexam`
  9. measure and subtract the median background in each image. Put down 100 random boxes around the image, find the *median* value and subtract. This takes advantage of the object masks created earlier to know where objects are.
  10. write out all the measured properties into a table so we know what happened and can see if anything went wrong.
8. Run `odi_scalestack_process.py` in the top-level folder. This will determine scaling factors, scale the image, and stack them into single images per filter.
  1. find sources in the image using `daofind`
  2. build a matched catalog of stars, making sure all the stars are visible on all the other images and not in a cell gap/etc.
  3. measure the fwhm to find an aperture size, then `iraf.apphot.phot` the sources in every image to get their total flux.
  4. choose a reference image (_LOWEST AIRMASS_)
  5. use the measured fluxes to compute a flux ratio between the reference image and all other images -- then divide each image by this ratio to scale each image to a common flux.
  6. stack the images using `iraf.imcombine(wcs='yes')`. This will take a long time since it is combining each individual OTA from each image. (with ODI5x6 data, around 270 images!) This produces a fits file named `object_filter.fits`.
9. Your images are now stacked. If you want to use our typical photometric calibration procedure, run `odi_phot_process.py` in the top-level folder. This will also (eventually) set some header keywords that will be useful to you in the future.
  1. Fix the WCS for the full image (one iteration). If things are good, this will only be a few arcseconds of shift.
  2. set the effective airmass for the stacked image to be that of the scaling reference image.
  3. get photometry of all SDSS catalog stars from the stacked images using a large aperture (4.5x fwhm).
  4. do matched pair photometric calibration using the photometry of the SDSS stars and the catalog values. See `full_calibrate.py:388` for full details.
  5. determine an aperture correction for each filter using the brightest SDSS stars.
  6. find *all* sources in the images using `daofind`, `iraf.apphot.phot` them, and determine their RA and Dec.
  7. based on RA and Dec, match the sources between the images, and calibrate the photometry using predetermined coefficients. Print calibrated magnitudes out to a file, and make a CMD of the final sources.
