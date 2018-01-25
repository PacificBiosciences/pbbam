
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
* `CMake`_ (3.0+)
* `Google Test`_
* `htslib`_ (1.4+)

For additional languages:

* `SWIG`_ (3.0.5+)

For building API documentation locally:

* `Doxygen`_

For maximal convenience, install htslib and google test in the same parent directory you plan to install pbbam.

.. _Boost: http://www.boost.org/
.. _clang: http://clang.llvm.org/
.. _CMake: https://cmake.org/
.. _Meson: http://mesonbuild.com/
.. _Ninja: https://ninja-build.org/ (only required when using Meson, optional for CMake)
.. _Doxygen: http://www.stack.nl/~dimitri/doxygen/
.. _gcc: https://gcc.gnu.org/
.. _Google Test: https://github.com/google/googletest
.. _htslib: https://github.com/samtools/htslib.git 
.. _SWIG: http://www.swig.org/

.. _getting_started-build:

Clone & Build
-------------

.. note::

   The following steps are for building the C++ library and command-line utilities. 
   If you are integrating pbbam into a C#, Python, or R project, take a look at the 
   instructions for :ref:`additional languages <swig_bindings>`.

The basic steps for obtaining pbbam and building it from source are as follows:

Build and install htslib, per the project's instructions (or on OSX "brew install htslib").

Clone
^^^^^

You should first clone the repository:

.. code-block:: console

   $ git clone https://github.com/PacificBiosciences/pbbam.git
   $ cd pbbam

Building with CMake
^^^^^^^^^^^^^^^^^^^

When building with CMake, create a separate build directory:

.. code-block:: console

   $ mkdir build
   $ cd build
   $ cmake ..
   $ make -j 4    # compiles using 4 threads

Output:

  * Library   : <pbbam_root>/lib
  * Headers   : <pbbam_root>/include
  * Utilities : <pbbam_root>/bin
 
You may need to set a few options on the cmake command, to point to dependencies' install locations. 
Common installation-related options include:

  * GTEST_SRC_DIR
  
Add these using the '-D' argument, like this:

.. code-block:: console

   $ cmake .. -DGTEST_SRC_DIR="path/to/googletest"
 
To run the test suite, run:

.. code-block:: console

   $ make test

To build a local copy of the (Doxygen-style) API documentation, run:

.. code-block:: console

   $ make doc
   
And then open <pbbam_root>/docs/html/index.html in your favorite browser.

.. _getting_started-integrate:

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

CMake-based projects
````````````````````

For CMake-based projects that will "ship with" or otherwise live alongside pbbam, you can 
use the approach described here.

Before defining your library or executable, add the following:

.. code-block:: cmake

   add_subdirectory(<path/to/pbbam> external/build/pbbam)

When it's time to run "make" this will ensure that pbbam will be built, inside your own project's 
build directory. After this point in the CMakeLists.txt file(s), a few variables will be available 
that can be used to setup your include paths and library linking targets:

.. code-block:: cmake

   include_directories( 
       ${PacBioBAM_INCLUDE_DIRS} 
       # other includes that your project needs
   )

   add_executable(foo)
   
   target_link_libraries(foo 
       ${PacBioBAM_LIBRARIES}
       # other libs that your project needs
   )

Non-CMake projects
``````````````````

If you're using something other than CMake for your project's build system, then you need to point 
it to pbbam's include directory & library, as well as those of its dependencies (primarily htslib).

If you built and installed pbbam using Meson, pkg-config files will be available to be consumed by
projects wishing to utilize pbbam. Autoconf, CMake, Waf, SCons and Meson all have means to determine
dependency information from pkg-config files.
