Running GPU-accelerated ``tblite`` calculations
================================================ 

If `tblite` has been compiled with GPU support (see :ref:`Installation <install>`), you can run GPU-accelerated calculations.
Currently, only the density building step is offloaded to the GPU, while other parts of the calculation (e.g., Hamiltonian construction, energy evaluation, etc.) are still performed on the CPU.

The `GAMBITS`_ library implements a number of alterntive density build algorithms, based on density matrix purification (DMP) techniques mainly by Niklasson and co-workers \ :footcite:`niklasson2002,niklasson2003,rubensson2014`.
These algorithms are particularly well suited for GPU acceleration, as they mainly involve matrix-matrix multiplications, which can be efficiently performed on GPUs.
In addition to these purification-based methods, a standard diagonalization-based density build method is also available based on NVIDIA's cuSOLVER library.

.. note::
   The implemented density matrix purification methods only support calculations with an electronic temperature of 0 K.
   If you need to perform calculations at finite electronic temperature, please use the diagonalization-based density build methods.
   A GPU-accelerated ``DSYGVD`` solver is also available, see below.

.. tip::
    The ``sp2`` algorithm is recommended for most applications due to its numerical stability and efficiency.
    By default, systems with less than 750 basis functions are computed on the CPU if a DMP methods have been chosen.
    For larger systems, the GPU is used automatically.
    When using the GPU, by default, the `mixed-precision` mode presented in \ :footcite:`steinbach2025` is used.

More details on the implemented methods and numerical accuracy considerations can be found in the original publications \ :footcite:`steinbach2025`.

The following table gives an overview on the available solvers, which can be requested using the ``--solver <solver_type>`` command line option.

=========== ================================================= ==========================
 Keyword(s)     Description                                     Reference
=========== ================================================= ==========================
 ``sp2``        2nd order spectral projection algorithm           :footcite:`niklasson2002`
 ``sp2-accel``   accelerated-version of SP2 (less stable)          :footcite:`rubensson2014`
 ``trs4``        trace-resetting 4th order algorithm               :footcite:`niklasson2003`
 ``gvd``         generalized eigenvalue equation solver (DSYGVD)    LAPACK `SYGVD`_
 ``gvd-gpu``     GPU-accelerated GVD solver using cuSOLVER          cuSOLVER `DSYGVD`_
 ``gvr``         eigenvalue equation solver (DSYGVR)                LAPACK `SYEVR`_
=========== ================================================= ==========================

.. _SYGVD: https://netlib.org/lapack/explore-html/d5/d0e/group__hegvd_ga1ed76f7825e99ff52abc815b4bf80e41.html
.. _SYEVR: https://netlib.org/lapack/explore-html/d1/d56/group__heevr_gaa2c43f6fba2353962ff981b759e4acab.html
.. _DSYGVD: https://docs.nvidia.com/cuda/cusolver/index.html#cusolverdn-t-sygvd
.. _GAMBITS: https://git.rwth-aachen.de/bannwarthlab/gambits