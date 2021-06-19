Usage in DFTB+
==============

The *tblite* project started as an effort to make the xTB Hamiltonian available in the `DFTB+`_ program package.\ :footcite:`hourahine2020`
DFTB+ supports a communication bridge to the *tblite* library when setting ``-DWITH_TBLITE=true`` at the CMake configuration step.

.. note::

   Support in DFTB+ will be presumably available with version 21.2 release in fall 2021.

.. _dftb+: https://github.com/dftbplus/dftbplus


Input structure
---------------

.. code-block:: bash

   Hamiltonian = xTB {
     method = "GFN1-xTB"
   }


Literature
----------

.. footbibliography::
