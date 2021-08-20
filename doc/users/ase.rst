ASE calculator
==============

The Python API of *tblite* natively support integration with the atomic simulation environment (`ASE`_).
By constructing a calculator most functionality of ASE is readily available.
For details on building the Python API checkout the :ref:`installation guide <python-build>`.

.. _ase: https://wiki.fysik.dtu.dk/ase/


Creating an ASE calculator
--------------------------

An ASE calculator can be constructed by using the *TBLite* class provided by the *tblite.ase* module.
For example to perform a single point calculation for a CO\ :sub:`2` crystal use

.. code-block:: python

   from tblite.ase import TBLite
   from ase.atoms import Atoms
   import numpy as np

   atoms = Atoms(
       symbols="C4O8",
       positions=np.array(
           [
               [0.9441259872, 0.9437851680, 0.9543505632],
               [3.7179966528, 0.9556570368, 3.7316862240],
               [3.7159517376, 3.7149292800, 0.9692330016],
               [0.9529872864, 3.7220864832, 3.7296981120],
               [1.6213905408, 1.6190616096, 1.6313879040],
               [0.2656685664, 0.2694175776, 0.2776540416],
               [4.3914553920, 1.6346256864, 3.0545920000],
               [3.0440834880, 0.2764611744, 4.4080419264],
               [4.3910577696, 3.0416409504, 0.2881058304],
               [3.0399936576, 4.3879335936, 1.6497353376],
               [0.2741322432, 4.4003734944, 3.0573754368],
               [1.6312174944, 3.0434586528, 4.4023048032],
           ]
       ),
       cell=np.array([5.68032, 5.68032, 5.68032]),
       pbc=np.array([True, True, True]),
   )

   atoms.calc = TBLite(method="GFN1-xTB")

   print(atoms.get_potential_energy())  # -1257.0943962462964 eV

The resulting calculator can be used like most ASE calculator, *e.g.* for optimizing geometries.
