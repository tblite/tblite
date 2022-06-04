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


.. doxygenfile:: tblite.h

.. doxygenfile:: tblite/version.h

.. doxygenfile:: tblite/macros.h

.. doxygenfile:: tblite/error.h

.. doxygenfile:: tblite/context.h

.. doxygenfile:: tblite/structure.h

.. doxygenfile:: tblite/calculator.h

.. doxygenfile:: tblite/container.h

.. doxygenfile:: tblite/result.h

.. doxygenfile:: tblite/param.h

.. doxygenfile:: tblite/table.h
