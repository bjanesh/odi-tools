# odi-tools

A suite of tools written in Pyraf, Astropy, Scipy, and Numpy to process individual QuickReduced images into single stacked images using a set of "best practices" for ODI data.

### Installation
simply clone this repository onto your local machine: `git clone https://github.iu.edu/wjanesh/odi-tools.git`

optionally add this folder to your `$PATH`

to run the scripts you'll need to install a number of dependencies:

`pip install numpy scipy astropy photutils pyraf tqdm matplotlib pandas`

It is possible to install these packages without root access by using the `--user` option:

`pip install --user package-name`

As noted on the [astropy website](http://astropy.readthedocs.org/en/stable/install.html), it might also be beneficial to use the `--no-deps` 
option when installing astropy to stop pip from automatically upgrading any of your previously installed packages, such as numpy.

`pip install --no-deps astropy`

### Usage
All you need to do to get started is download your QR-ed data from the ODI-PPA using the `wget` download command, then follow these steps.

1. move all individual .fz files into the top level folder: `mv calibrated/**/*.fz .`
2. unpack the compressed fits files using `funpack` [link](https://heasarc.gsfc.nasa.gov/fitsio/fpack/)
3. you need to rename your files to match the appropriate dither pointing identification. for example, QR files are named by the pattern `OBSID_OBJECT_FILTER.JOBID.fits`. The final digit of the OBSID `e.g. 20151008T195949.1` needs to match the number sequence of the dithers 1-9. Your data may not match this pattern due to restarted observations, multiple night observations, etc. Support will be added in the future for FITS header keywords providing this functionality without user interaction.
4. in `odi_process.py`, lines 15, 18, 20, indicate the proper filter names that you would like to reduce.
5. run `odi_process.py` in the folder containing the unpacked/renamed fits images. This will illumination correct the images, fix their WCS, reproject them to a common pixel scale, and perform background subtraction on them.
6. this will take a while, so make sure nothing bad happened
7. repeat step 4 by modifying the appropriate lines in `odi_scalestack_process.py`
8. run `odi_scalestack_process.py` in the folder containing the unpacked/renamed fits images. This will detect bright stellar sources in the images and use them to calculate a scaling factor relative to the image in the sequence with the lowest airmass, then apply the scale, stack the images, then add in a common background value.
9. finished! check your images to make sure everything went okay. ODI 5x6 data is currently marginally supported, and development is ongoing.
