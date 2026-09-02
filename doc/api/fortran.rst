Fortran API
===========

The *tblite* library seamlessly integrates with other Fortran projects via module interfaces,

.. note::

   Generally, all quantities used in the library are stored in `atomic units <https://en.wikipedia.org/wiki/Hartree_atomic_units>`_.

.. toctree::

   Full reference <https://tblite.github.io/tblite>


Handling of geometries and structure
------------------------------------

The basic infrastructure to handle molecular and periodic structures is provided by the `modular computation tool chain library <https://github.com/grimme-lab/mctc-lib>`_.
The library provides a structure type which is used to represent all geometry related informations in *tblite*.
A structure type can be constructed from arrays or read from a file.

The constructor is provided with the generic interface ``new`` and takes an array of atomic numbers (``integer``) or element symbols (``character(len=*)``) as well as the cartesian coordinates in Bohr.
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


Error handling
--------------

The basic error handler is an allocatable derived type, available from ``mctc_env`` as ``error_type``, which signals an error by its allocation status.

.. code-block:: fortran

   use mctc_env, only : error_type, fatal_error
   implicit none
   type(error_type), allocatable :: error

   call always_ok(error)
   if (allocated(error)) then
      print '(a)', "Unexpected failure:", error%message
   end if

   call always_failed(error)
   if (allocated(error)) then
      print '(a)', "Error:", error%message
   end if

   contains
      subroutine always_ok(error)
         type(error_type), allocatable, intent(out) :: error
      end subroutine always_ok

      subroutine always_failed(error)
         type(error_type), allocatable, intent(out) :: error

         call fatal_error(error, "Message associated with this error")
      end subroutine always_failed
   end

An unhandled error might get dropped by the next procedure call.


Calculation context
-------------------

The calculation context is available with the ``context_type`` from the ``tblite_context`` module.
The context stores error messages generated while running which can be queried using the type bound function ``failed``.
To access the actual errors the messages can be removed using the type bound subroutine ``get_error``.

An output verbosity is available in the context as the member verbosity, all procedures with access to the context will default to the verbosity of the context unless the verbosity level is overwritten by an argument.
To cutomize the output the ``context_logger`` abstract base class is available.
It must implement a type bound ``message`` procedure, which is used by the context to create output.
This type can be used to create callbacks for customizing or redirecting the output of the library.

The context also carries the work partition of the calculation, see :ref:`work-partition`.


.. _work-partition:

Work partitioning
-----------------

The ``tblite_partition`` module provides the ``work_partition`` type, which assigns a disjoint share of the interaction loops to each part of a distributed calculation.
Parts are zero based, every unit of work belongs to exactly one part, and summing the contributions of all parts reproduces the complete result.
*tblite* performs no communication itself, the reduction is left to the caller.

A partition is created with ``new_work_partition``, out of range parts are reported in the error handler.
The default constructed partition owns the complete work and is equivalent to not partitioning at all.

.. code-block:: fortran

   use mctc_env, only : error_type
   use tblite_partition, only : work_partition, new_work_partition
   implicit none
   type(error_type), allocatable :: error
   type(work_partition) :: partition

   ! every rank holds the complete structure and evaluates its own share
   call new_work_partition(error, partition, rank, nranks)

The partition is applied to a calculator with the type bound ``set_partition`` procedure, which propagates it to every interaction container, including containers added later with ``push_back``.
Alternatively the partition can be stored in the calculation context with ``ctx%set_partition(part, nparts, error)`` and handed to the calculator from there.

.. code-block:: fortran

   call calc%set_partition(partition)
   call xtb_singlepoint(ctx, mol, calc, wfn, accuracy, energy, gradient, sigma)

   ! tblite performs no communication, the caller reduces the partial results
   call mpi_allreduce(MPI_IN_PLACE, energy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm)
   call mpi_allreduce(MPI_IN_PLACE, gradient, size(gradient), MPI_DOUBLE_PRECISION, MPI_SUM, comm)

.. note::

   Structure dependent quantities such as coordination numbers, Born radii and the interaction caches are evaluated for the full system on every part.
   Only the interaction loops are partitioned, so the speedup is bound by those loops.

The diatomic blocks of the overlap, multipole and core Hamiltonian integrals and of the Hamiltonian gradient are partitioned as well, following the entries of the neighbour list rather than the atom pairs, so the share of each part is even for sparse and periodic systems.
The Born interaction matrix of the ALPB/GBSA model and the solvent accessible surface of the CDS term are partitioned over atom pairs and atoms, respectively.

Contributions which are not expressible as an interaction loop are carried in full by the first part.
This currently applies to the D3 and D4 dispersion corrections, which cannot partition their own loops yet, to the ddX solvation models, to the analytical linearized Poisson-Boltzmann gradient, whose inertia tensor couples all atoms, and to the external electric field.
The Born radii themselves enter non-linearly and are evaluated for the full system on every part.

Because the potential shifts of the self-consistent containers are partitioned as well, a partitioned calculation is only self-consistent if the potential is reduced in every iteration.
Either let *tblite* do this over MPI, see :ref:`mpi`, or use the partition on the individual containers and building blocks rather than on the full self-consistent driver.


.. _mpi:

Distributing over MPI
---------------------

MPI support is opt-in and has to be requested at build time with ``-Dmpi=true`` (meson) or ``-DTBLITE_WITH_MPI=ON`` (CMake).
Whether a build supports it can be queried at compile time with the ``tblite_has_mpi`` parameter and at runtime with ``get_tblite_feature("mpi")``, both from the ``tblite_features`` module.
Without MPI support every entry point of the ``tblite_mpi_utils`` module reports an error instead of performing communication.

.. code-block:: fortran

   use tblite_features, only : tblite_has_mpi, get_tblite_feature

   if (.not.get_tblite_feature("mpi")) error stop "tblite was built without MPI support"

With MPI enabled the calculation context can distribute the interaction loops over a communicator and reduce the partial results inside the library.
``set_mpi`` derives the work partition from the rank and size of the communicator, which defaults to ``MPI_COMM_WORLD``.
The partition still has to be handed to the calculator, a calculator that does not share the partition of the context is rejected rather than silently double counting or dropping contributions.

.. code-block:: fortran

   call mpi_init(stat)

   call ctx%set_mpi(error)                ! or ctx%set_mpi(error, comm)
   call calc%set_partition(ctx%partition)

   ! energy, gradient and virial are already reduced, every rank holds the total
   call xtb_singlepoint(ctx, mol, calc, wfn, accuracy, energy, gradient, sigma)

Internally *tblite* uses the ``mpi_f08`` interfaces, but communicators cross the library boundary as plain integer handles so that no MPI types leak into the calculation context or the calculator.
Users of ``mpi_f08`` pass ``comm%MPI_VAL``, users of the older ``mpi`` module pass the communicator directly.

The library reduces the density dependent potential in every self-consistent iteration, so all ranks follow the same SCF trajectory and end up with the same wavefunction.
The integral and core Hamiltonian matrices are reduced once after they are built, the diagonalization is then performed redundantly on every rank.
A failure on any rank is made visible to all of them, a rank leaving a collective on its own would deadlock the remaining ones.
``ceh_singlepoint`` supports the same distribution.
*tblite* neither initializes nor finalizes MPI, this remains the responsibility of the caller.

The ``tblite`` command line driver does this for you.
An MPI enabled binary always enters the MPI environment and derives its work partition from ``MPI_COMM_WORLD``, so running it under ``mpiexec`` distributes the calculation without any further option.
A single rank owns the complete work, which makes a normal invocation behave exactly as before.
Only the first rank reports and writes result files, every rank holds the same reduced result.

.. code-block:: shell

   mpiexec -n 4 tblite run --method gfn2 --grad struc.xyz

.. note::

   The xTB-ML features are rejected for a distributed calculation.
   They are evaluated from the partitioned interaction caches and are normalized by the total energy, so the partial results of the ranks cannot be summed afterwards.
   Bond orders and multipole moments are computed from the reduced wavefunction and remain available.

.. note::

   Reducing the integral matrices costs :math:`\mathcal{O}(N_\text{ao}^2)` communication per geometry and every rank still holds the full matrices, so memory does not scale with the number of ranks.


High-level interface
--------------------

The high-level interface is defined by the calculation context, the calculator instance and its restart data.
The calculation context is defined with the ``context_type``, which stores general settings regarding the overall method independent setup of the calculation.
The actual parametrisation data is stored in the ``xtb_calculator`` type.
An instance of the calculator can be used in a thread-safe way to perform calculations for a specific structure (defined by its number of atoms, unique elements and boundary conditions).
Changing the specific structure parameters requires to reconstruct the calculator.
Finally the specific persient data for a geometry is stored in a ``wavefunction_type``, which allows to restart calculations based on previous results.
