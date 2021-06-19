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

The interaction potential is parametrized by a Klopman–Ohno type potential in the xTB Hamiltonian.


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
