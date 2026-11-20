Running GPU-accelerated ``tblite`` calculations
===============================================

For a single-GPU GFN2-xTB calculation, enable the NVIDIA eigensolver at build time
(see :ref:`Installation <install>`) and select ``gvd-gpu`` on the command line.
For example, to compute an energy and gradient for an XYZ structure, run

.. code-block:: console

   $ tblite run --method gfn2 --solver gvd-gpu --grad gradient.txt structure.xyz

cuSOLVER diagonalizes the generalized eigenproblem and cuBLAS constructs the density
and energy-weighted density. Hamiltonian construction and the remaining GFN2-xTB
terms run on the CPU. This solver supports finite electronic temperatures and
retains eigenvectors for properties and restart files.

The optional `GAMBITS`_ library implements density matrix purification (DMP) algorithms
developed by Niklasson and co-workers :footcite:`niklasson2002,niklasson2003,rubensson2014`.
These algorithms are particularly well suited for GPU acceleration, as they mainly involve matrix-matrix multiplications, which can be efficiently performed on GPUs.

Purification requires an electronic temperature of 0 K. To select SP2 explicitly
on the GPU, build with CUDA and GAMBITS support and run

.. code-block:: console

   $ tblite run --method gfn2 --solver sp2 --etemp 0 --purification-runmode gpu structure.xyz

With automatic device selection, purification runs on the CPU below 750 basis functions
and on the GPU from 750 basis functions in a CUDA-enabled build.
The default GPU precision is mixed :footcite:`steinbach2025`.

.. list-table:: Solver selection with ``--solver``
   :header-rows: 1

   * - Keyword
     - Algorithm
   * - ``sp2``
     - Second-order spectral projection :footcite:`niklasson2002`
   * - ``sp2-accel``
     - Accelerated SP2 :footcite:`rubensson2014`
   * - ``trs4``
     - Fourth-order trace-resetting purification :footcite:`niklasson2003`
   * - ``gvd``
     - CPU divide-and-conquer diagonalization, LAPACK `SYGVD`_
   * - ``gvd-gpu``
     - Single-GPU divide-and-conquer diagonalization, cuSOLVER `DSYGVD`_
   * - ``gvr``
     - CPU diagonalization using LAPACK `SYEVR`_

.. footbibliography::

.. _SYGVD: https://netlib.org/lapack/explore-html/d5/d0e/group__hegvd_ga1ed76f7825e99ff52abc815b4bf80e41.html
.. _SYEVR: https://netlib.org/lapack/explore-html/d1/d56/group__heevr_gaa2c43f6fba2353962ff981b759e4acab.html
.. _DSYGVD: https://docs.nvidia.com/cuda/cusolver/index.html#cusolverdn-t-sygvd
.. _GAMBITS: https://git.rwth-aachen.de/bannwarthlab/gambits
