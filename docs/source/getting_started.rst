
.. _getting_started:

Getting Started
===============

.. _getting_started-requirements:

Requirements
------------

These components will almost certainly already be on your system. 
 
* `gcc`_ (4.8+) OR `clang`_ (v3.1+)
* pthreads
* zlib

Double-check your compiler version, to be sure it is compatible.

.. code-block:: console

   $ g++ -v    
   $ clang -v  

Additional requirements:

* `Boost`_ (1.55+)
* `Meson`_ (0.48+)
* `Google Test`_
* `htslib`_ (1.4+)

For building API documentation locally:

* `Doxygen`_

For maximal convenience, install htslib and google test in the same parent directory you plan to install pbbam.

.. _Boost: http://www.boost.org/
.. _clang: http://clang.llvm.org/
.. _Meson: https://mesonbuild.com
.. _Ninja: https://ninja-build.org/ (only required when using Meson, optional for CMake)
.. _Doxygen: http://www.stack.nl/~dimitri/doxygen/
.. _gcc: https://gcc.gnu.org/
.. _Google Test: https://github.com/google/googletest
.. _htslib: https://github.com/samtools/htslib.git 

.. _getting_started-build:

Clone & Build
-------------

.. note::

   The following steps are for building the C++ library and command-line utilities. 

The basic steps for obtaining pbbam and building it from source are as follows:

Build and install htslib, per the project's instructions (or on OSX "brew install htslib").

Clone
^^^^^

You should first clone the repository:

.. code-block:: console

   $ git clone https://github.com/PacificBiosciences/pbbam.git
   $ cd pbbam

Building with Meson
^^^^^^^^^^^^^^^^^^^

Building with Meson is generally faster and more versatile. Meson strictly requires building out of source:

.. code-block:: console

   $ mkdir build
   $ cd build
   $ meson --prefix /my/install/prefix -Dtests=true ..
   $ ninja

where ninja will by default utilize a number of threads for compilation equal to the number of logical
cores on your system. Here ``-Dtests=true`` enables pulling in dependencies for testing. In
order to run the test suite, run:

.. code-block:: console

   $ ninja test

If you wish to install pbbam, run:

.. code-block:: console

   $ ninja install

and ninja will install pbbam to ``/my/install/prefix``.

Integrate
---------

If you built and installed pbbam, pkg-config files will be available to be consumed by projects
wishing to utilize pbbam. Autoconf, CMake, Waf, SCons and Meson all have means to determine
dependency information from pkg-config files.
