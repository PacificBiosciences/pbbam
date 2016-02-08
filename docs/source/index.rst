.. pbbam documentation master file, created by
   sphinx-quickstart on Fri Dec  4 10:08:52 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _home:

pbbam documentation
===================

As of the 3.0 release of SMRTanalysis, PacBio is embracing the industry standard BAM 
format for (both aligned and unaligned) basecall data files. We have also formulated 
a BAM companion file format (bam.pbi) enabling fast access to a richer set of per-read 
information as well as compatibility for software built around the legacy cmp.h5 format.

The **pbbam** software package provides components to create, query, & edit PacBio BAM
files and associated indices. These components include a core C++ library, bindings for 
additional languages, and command-line utilities.

.. toctree::
   :maxdepth: 1

   getting_started
   api_reference
   swig_bindings
   commandline_utilities


Search:

* :ref:`genindex`
* :ref:`search`

