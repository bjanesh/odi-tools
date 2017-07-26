Measuring Stellar FWHM
======================

There are a number of steps in ``odi-tools`` that require having a measurement
of the stellar fwhm of sources on individual OTAs or on a fully stacked image.
In order to get these measurements we use the ``pyraf`` task ``rimexam`` on
a list of known x and y positions for SDSS sources on a given field.
Here is how the parameters are set for ``rimexam``::

    iraf.tv.rimexam.setParam('radius',radius)
    iraf.tv.rimexam.setParam('buffer',buff)
    iraf.tv.rimexam.setParam('width',width)
    iraf.tv.rimexam.setParam('rplot',20.)
    iraf.tv.rimexam.setParam('center','yes')
    iraf.tv.rimexam.setParam('fittype','gaussian')
    iraf.tv.rimexam.setParam('iterati',1)

.. automodule:: odi_helpers
   :members: getfwhm_ota, getfwhm_full
   :show-inheritance:
