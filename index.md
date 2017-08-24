---
---
### Installation
simply fork this repository and clone onto your local machine, e.g.: `git clone https://github.com/bjanesh/odi-tools.git`

optionally add this folder to your `$PATH`

to run the scripts you'll need to install a number of dependencies:

`pip install numpy scipy astropy astroquery pyraf matplotlib pandas photutils tqdm`

It is possible to install these packages without root access by using the `--user` option:

`pip install --user package-name`

As noted on the [astropy website](http://astropy.readthedocs.org/en/stable/install.html), it might also be beneficial to use the `--no-deps`
option when installing astropy to stop pip from automatically upgrading any of your previously installed packages, such as numpy.

`pip install --no-deps astropy`

### Basic Usage
All you need to do to get started is download your QR-ed data from the ODI-PPA using the `wget` download command, then follow these steps (if you aren't sure what to do, see below for more details).

1. rename the `images` folder created by the `wget` command to avoid confusion
2. unpack the compressed fits files using `funpack` [link](https://heasarc.gsfc.nasa.gov/fitsio/fpack/)
3. copy `example_config.yaml` to your data directory as `config.yaml` and edit the file to match your preferences/data. You do not need to rename your images or use any particular numbering sequence, though using dither sequence ID numbers is recommended for clarity. You may use more than 9 images, but be sure to give each image a unique ID! For any additional dither sequences we recommend using, e.g., 11-19, 21-29, etc. **NOTE: the current version of the `example_config.yaml` file contains the current recommended set of options for `odi-tools`.** 
4. run `odi_process.py` in the folder containing the unpacked fits images. This will (optionally) illumination correct the images, reproject them to a common pixel scale, and perform background subtraction on them.
5. this will take a while, so make sure nothing bad happened
6. run `odi_scalestack_process.py` in the folder containing the unpacked fits images. This will detect bright stellar sources in the images and use them to calculate a scaling factor relative to the image in the sequence with the lowest airmass, then apply the scale, stack the images, then add in a common background value. Finally, the images are flipped and optionally, aligned using integer pixel shifts.
7. Finished! Check your images to make sure everything went okay.
