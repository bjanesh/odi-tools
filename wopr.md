---
---
`wopr> Greetings, Professor Falken.`

`wopr> Would you like to play a game of ODI data processing?`

**Summary** 

`wopr` is a local IU Astro workstation. Its primary use is ODI data processing after the data has been reduced by the [QuickReduce pipeline](https://portal.odi.iu.edu). To account for a balance between processing speed and data volume, `wopr` is equipped with a 500GB solid state drive, a 1 TB spinning drive, 16GB of memory, and a fast workstation processor. We have also written a pipeline to process and stack pODI/5x6ODI data called [odi-tools](https://github.com/bjanesh/odi-tools), which is in use on `wopr`. Documentation for ODI data processing with `odi-tools` can be found [here](http://odi-tools.readthedocs.io).

**General Policies**

1. `wopr` **is a restricted access machine**, meaning you need a user account to use it. Contact Bob Lezotte to request an account. 

2. **Time on `wopr` is scheduled on a first-come-first-serve basis.** Schedule as much time as you need, but please respect others' time and use the scheduled time. Full processing of a 9-point dither pattern takes about 3-4 hours per filter, so you shouldn't need more than 2 days per target. **After your scheduled time, add one (1) day of "Bob" time** so that the computer can be reset for the next user. The schedule can be found at the bottom of this page or at [https://teamup.com/ks78bf366c93189e18](https://teamup.com/ks78bf366c93189e18). 

3. **All disk space on `wopr` will be treated as scratch space.** ODI data processing is extremely space intensive. The SSD is appropriately sized for full processing of ~30 5x6 ODI images. Use your space carefully! At the end of your scheduled time, **the SSD will be cleared.** Take your final data products with you! A few (representative) statistics about data size:
    * A single compressed 5x6 ODI image: `250 MB`
    * A single uncompressed 5x6 ODI image: `2.0 GB`
    * One 9-point dither pattern (unprocessed): `18 GB`
    * One 9-point dither pattern (processed): `130 GB`
    * Two 9-point dither patterns (processed): `260 GB`
    * Three 9-point dither patterns (processed): `390 GB`
    * Maximum SSD disk capacity: `490 GB`

4. **Final data products will be stored/archived on the IU Scholarly Data Archive.** The details of this process are still pending.

5. **Need help?** If you have technical issues with `wopr` please report them to Bob Lezotte. Critical issues with `odi-tools` should be reported to Bill Janesh or Owen Boberg. Non-critical issues or feature requests should be submitted as a [Github issue](https://github.com/bjanesh/odi-tools/issues). All other communication should be directed to the `wopr-l` mailing list.

**Getting started**

0. Log in to your account.
1. First time users should initialize `IRAF` by doing the following:
    1. navigate to your home directory: `cd ~`
    2. create an IRAF folder: `mkdir iraf`
    3. navigate to the iraf folder: `cd iraf`
    4. initialize IRAF: `mkiraf`
    5. (optional) copy your preferred `login.cl` file with `scp` or `sftp`
    6. start, and then exit pyraf: `pyraf`, then `.exit`
2. Get your QuickReduced data using the Download option on ODI-PPA. You probably want to use the `wget` option, so leave the "Tar exposure directories" box **unchecked**. Navigate to `/ssd1` and paste the command.
3. Your images will download to a folder named, appropriately, `images`. Rename this folder to match your object name, e.g. `mv images m15`.
4. Navigate into the folder with the images and uncompress them: `funpack *.fz`. This will take a few minutes. If you're worried about space, delete the compressed images `rm *.fz`. Otherwise you may wish to keep them until you are sure your data processed correctly.
5. Copy and edit a `config.yaml` file to include the appropriate options for your data: `cp $ODI_CONFIG config.yaml`. **UPDATE 26 January 2017:** illumination correction should be turned off.
6. Run the `odi-tools` data processing scripts `odi_process.py`, `odi_scalestack_process.py` as needed. 
7. Run the clean-up script `odi_cleanup.py`. This will *compress & move* your final products to a `.tar.gz` file in `/ssd1`. You can then copy this file to your own machine. Don't do this until you are satisfied with your final results!
8. Notify an administrator (Bob) when you are finished processing your data.

<iframe src="https://teamup.com/ks78bf366c93189e18" frameborder="0" width="100%" height="400"></iframe> 
