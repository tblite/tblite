Usage in DFTB+
==============

The *tblite* project started as an effort to make the xTB Hamiltonian available in the `DFTB+`_ program package.\ :footcite:`hourahine2020`
DFTB+ supports a communication bridge to the *tblite* library when setting ``-DWITH_TBLITE=true`` at the CMake configuration step.

.. note::

   Support in DFTB+ will be presumably available with version 21.2 release in fall 2021.

.. _dftb+: https://github.com/dftbplus/dftbplus


Input structure
---------------

Input to DFTB+ is provided by creating a file named ``dftb_in.hsd`` in the custom human-friendly structured data (HSD) format, which is inspired by XML.
The geometry input is provided in the *Geometry* group, either as *xyzFormat* for molecular geometries or as *vaspFormat* for periodic geometries.
To enable the xTB methods the *Hamiltonian* group is set to *xTB*.
In the *xTB* group the *Method* keyword is provided to select between GFN2-xTB, GFN1-xTB and IPEA1-xTB.

.. code-block:: bash
   :caption: Molecular GFN2-xTB calculation

   Geometry = xyzFormat {
   <<< "struc.xyz"
   }

   Hamiltonian = xTB {
     Method = "GFN2-xTB"
   }

To perform periodic calculations with the xTB Hamiltonian only the k-point sampling has to be added to the *xTB* group with *kPointsAndWeights*.
Using the *SuperCellFolding* provides Monkhorst–Pack k-point sampling.

.. code-block:: bash
   :caption: Periodic GFN1-xTB calculation

   Geometry = vaspFormat {
   <<< "POSCAR"
   }

   Hamiltonian = xTB {
     Method = "GFN1-xTB"
     kPointsAndWeights = SuperCellFolding {
       2   0   0
       0   2   0
       0   0   2
       0.5 0.5 0.5
     }
   }

Instead of providing a *Method* the xTB method can be initialized from a parameter file by providing its path in *ParameterFile*.

.. code-block:: bash

   # ...
   Hamiltonian = xTB {
     ParameterFile = "gfn2-xtb.toml"
     # ...
   }

Finally, to perform more than just single point calculations, the *Driver* group has to be provided.
Possible geometry optimizers are *ConjugateGradient* or *LBFGS*.
For periodic structures the lattice optimization can be enabled by setting *LatticeOpt* to *Yes*.

.. code-block:: bash

   Driver = ConjugateGradient {
     LatticeOpt = Yes
   }

For the full capabilities for DFTB+ check the `reference manual <https://dftbplus.org>`_.


Literature
----------

.. footbibliography::
