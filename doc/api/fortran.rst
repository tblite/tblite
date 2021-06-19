Fortran API
===========

The *tblite* library seamlessly integrates with other Fortran projects via module interfaces,

.. note::

   Generally, all quantities used in the library are stored in `atomic units <https://en.wikipedia.org/wiki/Hartree_atomic_units>`_.

.. contents::


Handling of geometries and structure
------------------------------------

The basic infrastructure to handle molecular and periodic structures is provided by the `modular computation tool chain library <https://github.com/grimme-lab/mctc-lib>`_.
The library provides a structure type which is used to represent all geometry related informations in *tblite*.
A structure type can be constructed from arrays or read from a file.

The array constructor is provided with the generic interface ``new`` and takes an array of atomic numbers (``integer``) or element symbols (``character(len=*)``) as well as the cartesian coordinates in Bohr.
Additionally, the molecular charge and the number of unpaired electrons can be provided the ``charge`` and ``uhf`` keyword, respectively.
To create a periodic structure the lattice parameters can be passed as 3 by 3 matrix with the ``lattice`` keyword.

An example for using the constructor is given here

.. code-block:: fortran

   subroutine example
      use mctc_env, only : wp
      use mctc_io, only : structure_type, new
      implicit none
      type(structure_type) :: mol
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: num(:)

      num = [6, 1, 1, 1, 1]
      xyz = reshape([ &
        &  0.00000000000000_wp, -0.00000000000000_wp,  0.00000000000000_wp, &
        & -1.19220800552211_wp,  1.19220800552211_wp,  1.19220800552211_wp, &
        &  1.19220800552211_wp, -1.19220800552211_wp,  1.19220800552211_wp, &
        & -1.19220800552211_wp, -1.19220800552211_wp, -1.19220800552211_wp, &
        &  1.19220800552211_wp,  1.19220800552211_wp, -1.19220800552211_wp],&
        & [3, size(num)])

      call new(mol, num, xyz, charge=0.0_wp, uhf=0)

      ! ...
   end subroutine example


To interact with common input file formats for structures the ``read_structure`` procedure is available.
The file type is inferred from the name of the file automatically or if a file type hint is provided directly from the enumerator of available file types.
The ``read_structure`` routine can also use an already opened unit, but in this case the file type hint is mandatory to select the correct format to read from.

.. code-block:: fortran

   subroutine example
      use mctc_env, only : error_type
      use mctc_io, only : structure_type, read_structure, file_type
      implicit none
      type(structure_type) :: mol
      type(error_type), allocatable :: error
      character(len=:), allocatable :: input

      input = "struc.xyz"

      call read_structure(mol, input, error, file_type%xyz)
      if (allocated(error)) then
         print '(a)', error%message
         stop 1
      end if

      ! ...
   end subroutine example


The structure type as well as the error type are using only allocatable members and can therefore be used without requiring explicit deconstruction.

Certain members of the structure type should be considered immutable, like the number of atoms (``nat``), the identifiers for unique atoms (``id``) and the boundary conditions (``periodic``).
To change those specific structure parameters the structure type and all dependent objects should be reconstructed to ensure a consistent setup.
Other properties, like the geometry (``xyz``), molecular charge (``charge``), number of unpaired electrons (``uhf``) and lattice parameters (``lattice``) can be changed without requiring to reconstruct dependent objects like calculators or restart data.


High-level interface
--------------------

The high-level interface is defined by the calculation context, the calculator instance and its restart data.
The calculation context is defined with the ``context_type``, which stores general settings regarding the overall method independent setup of the calculation.
The actual parametrisation data is stored in the ``xtb_calculator`` type.
An instance of the calculator can be used in a thread-safe way to perform calculations for a specific structure (defined by its number of atoms, unique elements and boundary conditions).
Changing the specific structure parameters requires to reconstruct the calculator.
Finally the specific persient data for a geometry is stored in a ``wavefunction_type``, which allows to restart calculations based on previous results.
