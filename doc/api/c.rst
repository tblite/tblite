C API
=====

The C API bindings are provided by using the ``iso_c_binding`` intrinsic module.
Generally, objects are exported as opaque pointers and can only be manipulated within the library.
The API user is required delete all objects created in the library by using the provided deconstructor functions to avoid mamory leaks.

Overall four classes of objects are provided by the library

- error handlers (``tblite_error`` and ``tblite_context``),
  used to communicate exceptional conditions and errors from the library to the user,
  also usable to communicate general setting to the library in case of a context handle
- structure containers (``tblite_structure``),
  used to represent the system specific information and geometry data,
  only the latter are mutable for the user
- calculator objects (``tblite_calculator``),
  used to store actual method parametrisation
- result and restart data (``tblite_result``)
  used to hold calculated properties and persistent data for restarting between calculations

.. note::

   Generally, all quantities provided to the library are assumed to be in `atomic units <https://en.wikipedia.org/wiki/Hartree_atomic_units>`_.

.. contents::


Error handling
--------------

The library provides two different kinds handlers, a light error handle type (``tblite_error``) and a environment context type (``tblite_context``).
The error handle is used in the context of simple tasks and requires only small overhead to construct and use.
It is mainly used in the context of retrieving data or building structures.
The environment context can be considered a persistent setup for all calculations performed with the library, it is usually used together with calculator objects.
While the error handle can only contain a single error, multiple errors can be accumulated in a context object, which allows storing more complex error information like they can occur in an actual calculation.

Both handleres are represented by an opaque pointer and can only be manipulated by call from the library.
The user of those objects is required to delete the handlers again using the library provided constructors to avoid memory leaks.


Structure data
--------------

The structure data is used to represent the system of interest in the library.
It contains immutable system specific information like the number of atoms, the unique atom groups and the boundary conditions as well as mutable geometry data like cartesian coordinates and lattice parameters.
