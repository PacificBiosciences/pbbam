.. _swig_bindings:

Additional Languages
====================

pbbam uses SWIG to generate bindings for other languages. Currently this includes support for C#, Python, and R.

These bindings are disabled by default. See the entry below for your target language to configure pbbam & integrate
the bindings into your project.

.. _swig_bindings-csharp:

C#
------

Building
````````

To build the support for C#, you need to tell CMake to enable it before building:

.. code-block:: console

   $ cmake .. -DPacBioBAM_wrap_csharp
   $ make

The 'make' step will build relevant libraries/wrappers, and then run a simple program using them, 
as a quick sanity-check. 

After building, the libraries and wrappers can be found under the pbbam/lib/csharp directory. 

API Example
```````````

.. code-block:: c#

   using PacBio.BAM;

   namespace TestStuff
   {
       public class TestPbbam
       {
           public static void TestZmwQuery()
           {
               var d = new DataSet("foo.bam");
               var q = new ZmwQuery(new IntList {1, 2, 3}, d);
               var q2 = new ZmwQuery(new IntList { 14743 }, d);
               if (0 != q.Count() || 4 != q2.Count())
               {
                   throw new Exception("ZmwQuery not working");
               }
               Console.WriteLine("TestZmwQuery - OK!");
           }
       }
   }

.. _swig_bindings-python:

Python
------

Building
````````

To build the support for Python, you need to tell CMake to enable it:

.. code-block:: console

   $ cmake .. -DPacBioBAM_wrap_python
   $ make

The 'make' step will build relevant libraries/wrappers, and then run a simple program using them, 
as a quick sanity-check. 

After building, the libraries and wrappers can be found in the pbbam/lib/python directory. 
'make test' will also include some Python-side unit tests as well.

To use the PacBioBam module, you can set your PYTHONPATH before invoking your script:

.. code-block:: console

   $ PYTHONPATH="path/to/pbbam/lib/python" python myScript.py

Or otherwise configure your environment to find the PacBioBam module. 

API Example
```````````

.. code-block:: python

   import PacBioBam
   
   try:
       file = PacBioBam.BamFile('foo.bam')
       writer = PacBioBam.BamWriter('new.bam', file.Header())
       dataset = PacBioBam.DataSet(file)
       entireFile = PacBioBam.EntireFileQuery(dataset)
       for record in PacBioBam.Iterate(entireFile):
           writer.Write(record)
   except RuntimeError:
       # found error
   
Python-Specific Notes
`````````````````````
   
Iteration
.........

Iteration over dataset queries in Python will likely need to use the PacBioBam.Iterate() method. Thus
file iteration loops will look something like the following:

.. code-block:: python
       
   entireFile = PacBioBam.EntireFileQuery("input.bam")
   for record in PacBioBam.Iterate(entireFile):
       foo.bar(record)

Exception Handling
..................
   
Exceptions are used widely by the C++ library. To handle them from Python, you can use try blocks, looking for
any RuntimeError:

.. code-block:: python

   try:
       file = PacBioBam.BamFile("does_not_exist.bam")
   except RuntimeError: 
       print("caught expected error")
   
.. _swig_bindings-r:

R
------

Building
````````

To build the support for R, you need to tell CMake to enable it:

.. code-block:: console

   $ cmake .. -DPacBioBAM_wrap_r
   $ make
   
The 'make' step will build relevant libraries/wrappers, and then run a simple program using them, 
as a quick sanity-check. 

After building, the libraries and wrappers can be found in the pbbam/lib/R directory. 
'make test' will also include some R-side unit tests as well.   

To use the PacBioBam module in your script, nothing should be needed up front - simply invoke 'R' as normal. 
You'll do the dynamic load of the R module near the beginning of your script:

.. code-block:: r

   # load pbbam R library
   lib_path <- "path/to/pbbam/lib/R"
   pbbam_libname <- paste(lib_path, "PacBioBam",   sep="/")
   pbbam_wrapper <- paste(lib_path, "PacBioBam.R", sep="/")
   dyn.load(paste(pbbam_libname, .Platform$dynlib.ext, sep=""))
   source(pbbam_wrapper)
   cacheMetaData(1) 


API Example
```````````

.. code-block:: r

   # load pbbam R library
   lib_path <- "path/to/pbbam/lib/R"
   pbbam_libname <- paste(lib_path, "PacBioBam",   sep="/")
   pbbam_wrapper <- paste(lib_path, "PacBioBam.R", sep="/")
   dyn.load(paste(pbbam_libname, .Platform$dynlib.ext, sep=""))
   source(pbbam_wrapper)
   cacheMetaData(1)    
  
   # sample method
   copyFileAndFetchRecordNames <-function(inputFn, outputFn) {
	
       result <- tryCatch(
       {
           file   <- BamFile(inputFn)
           writer <- BamWriter(outputFn, file$Header())
           ds     <- DataSet(file)
            
           entireFile <- EntireFileQuery(ds)
           iter <- entireFile$begin()
           end  <- entireFile$end()
   			
           while ( iter$'__ne__'(end) ) {
               record <- iter$value()
                
               names_in <- c(names_in, record$FullName())
               writer$Write(record)
               iter$incr()
            }
            writer$TryFlush()
            return(names_in)
        },
        error = function(e) {
            # handle error 
            return(list())
        })
        return(result)
   }

R-Specific Notes
````````````````

Iteration
.........

To compare iterators, you'll need to explicitly use the '__eq__' or '__ne__' methods. Thus iterating over
a data query, will look something like this:

.. code-block:: r

   iter <- query$begin()
   end  <- query$end()
   while ( iter$'__ne__'(end) ) {
       record <- iter$value() 
       
       # do stuff with record
   }
   
operator[]
..........  
   
In C++, operator[] can be used in some classes to directly access elements in a sequence, e.g. Cigar string

.. code-block:: cpp

   CigarOperation op = cigar[0]; 
   
For the R wrapper, if you want to do the same sort of thing, you'll need to use the '__getitem__' method. 
Please note that these are **0-based** indices, not 1-based as in much of R. 

.. code-block:: r

   op <- cigar$'__getitem__'(0) 
   
Exception Handling
..................

Exceptions are used widely by the C++ library. To handle them from R, you can use the 'tryCatch' block, listening for 
'error' type exceptions.

 .. code-block:: r
 
    result <- tryCatch(
    {
        f <- BamFile("does_not_exist.bam") # this statement will throw
    },
    error = function(e) {
        print(paste("caught expected erorr: ",e))
    })
