
QuickReduce, ODI-PPA, and odi-tools
===================================

QuickReduce is a set of pure python tools to reduce data from ODI. QuickReduce
was created by Ralf Kotulla (UW Milwaukee, UW Madison) for the WIYN Observatory.
The source code for QuickReduce is now on `github
<https://github.com/WIYN-ODI/QuickReduce>`_. Documenation for the QuickReduce
pipeline is available at this `link.
<http://members.galev.org/rkotulla/research/podi-pipeline/>`_

The `ODI-PPA <https://portal.odi.iu.edu/index/front>`_ is the online portal
used to access, sort, and run QuickReduce on your ODI data.
Information for gaining access and using the portal can be found at these
`help pages.
<https://help.odi.iu.edu/display/
help/ODI+Pipeline%2C+Portal+and+Archive+-+Quick+Guide>`_

**odi-tools** is designed to work on images that have been processed using
QuickReduce and downloaded from the ODI-PPA.

Running QuickReduce
*******************

After you become familiar with the layout and operation of the ODI-PPA, you can
use these following steps to process your data. The options we list here are
what we have determined to be the best practices for running QuickReduce. It is
important to remember that these options can be data dependent. We will update
these options should the change.

* Add the images you wish to reduce to a collection in the ODI-PPA and from
  the collection action menu, choose QuickReduce. See the PPA help pages for
  information about creating collections.


* Run QuickReduce with the following options. The ``[X]`` denotes that the
  option is selected.

   * ``[X]`` WCS
   * ``[X]`` Photometry
   * ``[X]`` Fringe (i- and z-band only)
   * ``[X]`` Persistency
   * ``[X]`` Nonlinearity
   * ``[X]`` Pupil Ghost (calibrations)
   * ``[ ]`` Pupil Ghost (science)
   * ``[X]`` Use Bad Pixel Masks
   * ``[X]`` Cosmic Ray Removal, 3 iteration(s)

* You should give this job a meaningful name and check the ``Email me when
  when finished`` box.

* When you receive the email letting you know your job is complete, download
  the QuickReduce data to your hard drive using the ``Download Results`` button
  on the QuickReduce job page. You don't need to select any of the optional boxes,
  just name the job and click submit. Eventually a ``wget`` command will pop up.
  Copy it to your clipboard, navigate to the folder you want the data to go into,
  then paste the wget command in your command line.

* At this point you are ready to start running ``odi-tools``. See the
  :ref:`basic_usage` documentation for information on starting this
  process. See
  :ref:`qr_nav` for a brief tutorial on getting to know a QuickReduced image.
