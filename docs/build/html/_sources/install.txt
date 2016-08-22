
.. _install:

.. ====================
.. Installing odi-tools
.. ====================

Installation
------------

To use this code simply fork this repository and clone it onto your
local machine::

    $ git clone https://github.com/bjanesh/odi-tools.git
    $ cd odi-tools

Optionally add this folder to your ``$PATH`` so the odi-scripts maybe
used in any current working directory.

To run the scripts you will  need to install a number of dependencies::

    $ pip install numpy scipy astropy photutils pyraf tqdm matplotlib pandas

It is possible to install these packages without root access by using the
``--user`` option::

    $ pip install --user package-name

As noted on the `astropy website <http://astropy.readthedocs.org
/en/stable/install.html>`_, it might also be beneficial to use the ``--no-deps``
option when installing astropy to stop ``pip`` from automatically upgrading any
of your previously installed packages, such as ``numpy``::

    $ pip install --no-deps astropy
