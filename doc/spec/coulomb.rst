.. _coulomb:

Electrostatic interactions
==========================

The Coulomb energy in tight binding is described using density fluctions δρ.
The form of the electrostatic energy is defined by the scheme used to model the density fluctuations δρ.
For extended tight binding a multipole expansion is used.

.. math::

   \delta\rho = \sum_\kappa q_\kappa + \mu_\kappa + \theta_\kappa

where q:sub:`κ` is the orbital partial charge, μ:sub:`κ` the orbital dipole moment, and θ:sub:`κ` the (traceless) orbital quadrupole moment for the orbital *κ*.
The partial charges, dipole moments and quadrupole moments are computed from density matrix *P* by Mulliken population analysis.

.. math::

   q_\kappa = n_\kappa - \sum_\lambda S_\kappa\lambda P_\lambda\kappa,

where *S* is the overlap matrix and *n*:sub:`κ` the occupation number of the orbital in the atomic reference.
Practically, the partial charges for a shell *l* or atom *A* are used instead of the orbital partial charge:

.. math::

   q_l = \sum_{\kappa\in l} q_\kappa \quad\text{or}\quad q_\text{A} \sum_{\kappa\in\text{A}} q_\kappa

Using the shell *l* instead of the orbital *κ* is preferred due the spherical symmetric atomic reference.
For example computed chemical hardness values for the orbitals will be identical within a shell due to using spherical symmetric atomic densities.

Similarly the atomic dipole moment is computed by population analysis

.. math::

   \mu_\text{A} = \sum_{\kappa\in\text{A}}\sum_\lambda D_{\kappa\lambda,\text{A}} P_\lambda\kappa

Notably, *D* is the dipole moment integral, with the aufpunkt of the dipole operator on the position of atom A.
For the atomic quadrupole moment the Mulliken population is done in a similar way

.. math::

   \theta_\text{A} = \sum_{\kappa\in\text{A}}\sum_\lambda Q_\kappa\lambda^{\kappa} P_\lambda\kappa

Second order
------------

The electrostatic energy in extended tight binding as second order in the density fluctuations and truncated based on the distance dependence of the Coulombic interactions.
In GFN1-xTB only the leading term, *R*:sup:`-1` is considered while in GFN2-xTB the terms up to *R*:sup:`-3` are included.
For GFN1-xTB subsequently only isotropic electrostatic interactions between partial charges are included while for GFN2-xTB anisotropic terms up to dipole-dipole and charge-quadrupole are considered.

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

For the molecular case the interaction kernel for the dipole-charge interactions is expressed as

.. math::

   \gamma^{k}_\text{AB} =
   \left(1 + 6 * \left(\frac{f_\text{A} + f_\text{B}}{2R_\text{AB}}\right)^{d_3}\right)^{-1}
   \frac{R_{\text{AB},k}}{R_\text{AB}^3}

with *f*:sub:`A/B` being a critical radius to remove short-range multipole interactions and *d*:sub:`3` being the exponent for the damping function.
For the dipole-dipole and charge-quadrupole interactions the kernel is written as

.. math::

   \gamma^{kk'}_\text{AB} =
   \left(1 + 6 * \left(\frac{f_\text{A} + f_\text{B}}{2R_\text{AB}}\right)^{d_5}\right)^{-1}
   \frac{\delta_{kk'}R_{\text{AB}}^2-3R_{\text{AB},k}R_{\text{AB},k'}}{R_\text{AB}^5}

here a different exponent *d*:sub:`5` is used for damping the dipole-dipole and charge-quadrupole interactions at short-range.

Third order
-----------

The isotropic third-order contributions are included as the trace of the on-site shell-resolved Hubbard derivatives.

.. math::

   E_\text{IXC} =
   \frac13 \sum_\text{A} \sum_{l}
   \Gamma^l_\text{A} (q^l_\text{A})^3
