Running ``tblite`` in parallel
==============================

The ``tblite`` program uses shared memory OpenMP parallelisation, to calculate larger systems
an appropriate OMP stacksize must be provided, chose a reasonable large number by

.. code:: bash
  > export OMP_STACKSIZE=4G
  
.. note::

   Note that the memory requirement will increase with the system size *and* the number
   of requested threads.

To distribute the number of threads reasonable in the OpenMP section
it is recommended to use

.. code:: bash
  > export OMP_NUM_THREADS=<ncores>,1
You might want to deactivate nested OMP constructs by

.. code:: bash
  > export OMP_MAX_ACTIVE_LEVELS=1
.. tip::

   Most OpenMP regions allow to customize the scheduling by setting ``OMP_SCHEDULE``,
   for many threads the ``dynamic`` schedule has proven to give a good load-balance
   between all threads.

Depending on the linear algebra backend used when compiling ``tblite``, OpenMP threaded versions are available.
Usually, those backends repect the settings entered by ``OMP_NUM_THREADS``.
However, you can still adjust the parallelisation indivdually for the linear algebra backend.
For intel's Math Kernel Library, the environment variable is ``MKL_NUM_THREADS``.
For the OpenBLAS backend, use ``OPENBLAS_NUM_THREADS`` instead.
It is then exported for the current session as follows:

.. code:: bash
  > export MKL_NUM_THREADS=<ncores>
or respectively:

.. code:: bash
  > export OPENBLAS_NUM_THREADS=<ncores>
When computing large systems the limit of memory to be used for variables 
saved on the stack should be adjusted as it can lead to segmentation faults.
This can be achieved on UNIX systems (Linux and Mac) through the ``ulimit`` program, as so:

.. code:: bash
  > ulimit -s unlimited
Parallelisation using the python API
-------------------------------------

When using ``tblite``'s python API, the parallelisation behavior is also controlled via the aforementioned environment variables.
These variables can be set in the terminal before launching the python code containing the tblite calculations.
Another possibility is to set the varaibles from within the python code.
This can be achieved by the ``os.environ`` object, for details consider their `documentation <https://docs.python.org/3/library/os.html#os.environ>` for details.

To set up OpenMP in an analogue fashion as above:

.. code:: python
   import os
   os.environ['OMP_STACKSIZE'] = '3G'
   os.environ['OMP_NUM_THREADS'] = '<ncores>,1'
   os.environ['OMP_MAX_ACTIVE_LEVELS'] = '1'
The maximal stack size can also set from within python, we tested this using the `resource <https://docs.python.org/3/library/resource.html#resource-limits>` module.

To set the stack size to unlimited the following code snippet can be used:

.. code:: python
   import resource
   resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))