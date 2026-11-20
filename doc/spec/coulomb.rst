.. _coulomb:

Electrostatic interactions
==========================

.. contents::


Second order
------------


Isotropic electrostatics
~~~~~~~~~~~~~~~~~~~~~~~~

The isotropic electrostatic in a shell-resolved formulation is given by the parametrized Coulomb interaction between shellwise partial charges

.. math::

   E_\text{IES} =
   \frac12 \sum_{\text{A},\text{B}} \sum_{l,l'}^{s,p,d}
   q^{l}_\text{A} \gamma^{ll'}_\text{AB} q^{l'}_\text{B}

The interaction potential is parametrized by a Klopman–Ohno type potential in the xTB Hamiltonian or the γ-functional as used in the DFTB Hamiltonian.

Klopman–Ohno kernel
^^^^^^^^^^^^^^^^^^^

The interaction kernel for the Klopman–Ohno electrostatic is given by

.. math::

   \gamma^{ll'}_\text{AB} =
   \left(
   R_\text{AB}^g + f_\text{av}(\eta_A^l, \eta_B^{l'})^{-g}
   \right)^{-\frac1g}

where η:sub:`A/B` are the chemical hardness parameters of the respective shells and *g* is the exponent to manipulate the potential shape.

For three-dimensional periodic systems with :math:`g=2`, the kernel is evaluated using a generalized Ewald partition.\ :footcite:`buccheri2025`
For a lattice translation :math:`\mathbf T`, the Klopman--Ohno kernel has the binomial expansion

.. math::

   \left(\lvert\mathbf R_{AB}+\mathbf T\rvert^2
   +\eta_{Al,Bl'}^{-2}\right)^{-1/2}
   = \sum_{j=0}^{\infty}\binom{-1/2}{j}
   \frac{\eta_{Al,Bl'}^{-2j}}
   {\lvert\mathbf R_{AB}+\mathbf T\rvert^{2j+1}}.

Retaining the two long-range terms gives the exact partition

.. math::

   \begin{split}
   \gamma_{Al,Bl'}^{\text{PBC}}
   ={}& \sum_{\mathbf T}^{\prime}
   \left[
   \left(r_{AB,\mathbf T}^2+\eta_{Al,Bl'}^{-2}\right)^{-1/2}
   -r_{AB,\mathbf T}^{-1}
   +\frac12\eta_{Al,Bl'}^{-2}r_{AB,\mathbf T}^{-3}
   \right] \\
   &+S_1(\mathbf R_{AB})
   -\frac12\eta_{Al,Bl'}^{-2}S_3(\mathbf R_{AB}),
   \end{split}

where :math:`r_{AB,\mathbf T}=\lvert\mathbf R_{AB}+\mathbf T\rvert` and the prime excludes :math:`r_{AB,\mathbf T}=0`.
The residual in square brackets decays as :math:`r^{-5}` and is summed in real space.
With Ewald parameter :math:`\alpha=\sqrt{\pi}K`, the Coulomb lattice sum is

.. math::

   \begin{split}
   S_1(\mathbf R) ={}&
   \sum_{\mathbf T}^{\prime}
   \frac{\operatorname{erfc}(\alpha r_{\mathbf T})}{r_{\mathbf T}}
   +\frac{4\pi}{V}\sum_{\mathbf G\ne0}
   \frac{\exp[-G^2/(4\alpha^2)]}{G^2}
   \cos(\mathbf G\cdot\mathbf R) \\
   &-\delta_{\mathbf R,0}\frac{2\alpha}{\sqrt\pi},
   \end{split}

and the cubic lattice sum is

.. math::

   \begin{split}
   S_3(\mathbf R) ={}&
   \sum_{\mathbf T}^{\prime}\left[
   \frac{\operatorname{erfc}(\alpha r_{\mathbf T})}{r_{\mathbf T}^3}
   +\frac{2\alpha\exp(-\alpha^2r_{\mathbf T}^2)}
   {\sqrt\pi r_{\mathbf T}^2}\right] \\
   &+\frac{2\pi}{V}\sum_{\mathbf G\ne0}
   E_1\!\left(\frac{G^2}{4\alpha^2}\right)
   \cos(\mathbf G\cdot\mathbf R) \\
   &+\frac{4\pi}{V}\left[
   \ln\!\left(\frac{\alpha}{\sqrt\pi}\right)
   +\frac12\left(\ln\pi-\psi\!\left(\frac32\right)\right)
   \right]
   -\delta_{\mathbf R,0}\frac{4\pi}{3}
   \left(\frac{\alpha}{\sqrt\pi}\right)^3.
   \end{split}

Here :math:`V` is the unit-cell volume, :math:`\mathbf G` is a reciprocal lattice vector, :math:`E_1(x)=\Gamma(0,x)` is the exponential integral, and :math:`\psi` is the digamma function.
The third line of :math:`S_3` contains the :math:`\mathbf G=0` contribution; it is required to make the result independent of :math:`\alpha`.
The final terms in :math:`S_1` and :math:`S_3` remove the Gaussian self-interaction.


γ-functional kernel
^^^^^^^^^^^^^^^^^^^

The interaction kernel for the DFTB γ-functional is derived from the integral of two exponential densities

.. math::

   \begin{split}
   \gamma^{ll'}_\text{AB} =
   \frac1{R_\text{AB}}
   - \exp[-\tau_\text{A}R]
     \left(
     \frac{\tau_\text{B}^4\tau_\text{A}}{2(\tau_\text{A}^2-\tau_\text{B}^2)^2}
     - \frac{\tau_\text{B}^6\tau_\text{A} - 3\tau_\text{B}^4\tau_\text{A}^2}
       {(\tau_\text{A}^2-\tau_\text{B}^2)^3 R_\text{AB}}
     \right)
     \\
   - \exp[-\tau_\text{B}R]
     \left(
     \frac{\tau_\text{A}^4\tau_\text{B}}{2(\tau_\text{B}^2-\tau_\text{A}^2)^2}
     - \frac{\tau_\text{A}^6\tau_\text{B} - 3\tau_\text{A}^4\tau_\text{B}^2}
       {(\tau_\text{B}^2-\tau_\text{A}^2)^3 R_\text{AB}}
     \right)
   \end{split}

where τ:sub:`A/B` are scaled Hubbard parameters of the respective shells and *R* is the distance between the atomic sides.


Anisotropic electrostatics
~~~~~~~~~~~~~~~~~~~~~~~~~~

The anisotropic electrostatic in an atom-resolved formulation is given by the multipole interactions between the different moments:

.. math::

   E_\text{AES} =
   \sum_{\text{A},\text{B}} \sum_{k}^{x,y,z}
   q_\text{A} \gamma^{k}_\text{AB} \mu^{k}_\text{B}
   + \frac12 \sum_{\text{A},\text{B}} \sum_{k,k'}^{x,y,z}
   \mu^{k}_\text{A} \gamma^{kk'}_\text{AB} \mu^{k'}_\text{B}
   + \sum_{\text{A},\text{B}} \sum_{k,k'}^{x,y,z}
   q_\text{A} \gamma^{kk'}_\text{AB} \theta^{kk'}_\text{B}


Third order
-----------

The isotropic third-order contributions are included as the trace of the on-site shell-resolved Hubbard derivatives.

.. math::

   E_\text{IXC} =
   \frac13 \sum_\text{A} \sum_{l}
   \Gamma^l_\text{A} (q^l_\text{A})^3


Literature
----------

.. footbibliography::