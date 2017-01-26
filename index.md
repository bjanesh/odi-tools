---
---
### Installation
simply fork this repository and clone onto your local machine, e.g.: `git clone https://github.com/bjanesh/odi-tools.git`

optionally add this folder to your `$PATH`

to run the scripts you'll need to install a number of dependencies:

`pip install numpy scipy astropy photutils pyraf tqdm matplotlib pandas`

It is possible to install these packages without root access by using the `--user` option:

`pip install --user package-name`

As noted on the [astropy website](http://astropy.readthedocs.org/en/stable/install.html), it might also be beneficial to use the `--no-deps`
option when installing astropy to stop pip from automatically upgrading any of your previously installed packages, such as numpy.

`pip install --no-deps astropy`

### Usage
All you need to do to get started is download your QR-ed data from the ODI-PPA using the `wget` download command, then follow these steps (if you aren't sure what to do, see below for more details).

1. move all individual .fz files into the top level folder: `mv calibrated/**/*.fz .`
2. unpack the compressed fits files using `funpack` [link](https://heasarc.gsfc.nasa.gov/fitsio/fpack/)
3. you need to rename your files to match the appropriate dither pointing identification. for example, QR files are named by the pattern `OBSID_OBJECT_FILTER.JOBID.fits`. The final digit of the OBSID `e.g. 20151008T195949.1` needs to match the number sequence of the dithers 1-9. Your data may not match this pattern due to restarted observations, multiple night observations, etc.
4. copy `example_config.yaml` to your data directory as `config.yaml` and edit the file to match your preferences/data. Make sure that the number for each image matches the correct number in the dither sequence!
5. run `odi_process.py` in the folder containing the unpacked/renamed fits images. This will (optionally) illumination correct the images, fix their WCS, reproject them to a common pixel scale, and perform background subtraction on them.
6. this will take a while, so make sure nothing bad happened
7. run `odi_scalestack_process.py` in the folder containing the unpacked/renamed fits images. This will detect bright stellar sources in the images and use them to calculate a scaling factor relative to the image in the sequence with the lowest airmass, then apply the scale, stack the images, then add in a common background value.
8. finished! check your images to make sure everything went okay.
