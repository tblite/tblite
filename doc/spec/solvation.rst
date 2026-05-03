.. _solvation:


Solvation
=========

Implicit solvation models are available for the calculation of the solvation free energies partitioned as 

.. math::
   \Delta G_{\text{solv}} = \Delta G_{\text{polar}} + \Delta G_{\text{npol}} + \Delta G_{\text{shift}}

including the polar contribution :math:`{\Delta G_{\text{polar}}}` based on electrostatics, the non-polar contribution :math:`{\Delta G_{\text{npol}}}` based on cavity formation and dispersion, and a constant shift :math:`{\Delta G_{\text{shift}}}` depending on the thermodynamic state of initial gas and final liquid solution.

ALPB solvation model
--------------------

The analytical linearized Poisson-Boltzmann (ALPB) model evaluates the polar contribution

.. math::
   \Delta G^{\text{ALPB}}_{\text{polar}} = 
   - \frac{1}{2} \left(\frac{1}{\epsilon_{\text{in}}} - \frac{1}{\epsilon_{\text{out}}}\right) 
   \frac{1}{1+\alpha\beta}
   \sum_{A,B} q_{A} q_{B} \left( \frac{1}{f(R_{AB, a_{A}, a_{B}})} + \frac{\alpha\beta}{\mathcal{A}_{\text{det}}} \right)

based on the ALPB constant :math:`{\alpha}` (set to 0.571214), the solute (:math:`{\epsilon_{\text{in}}=1}`) and solvent (:math:`{\epsilon_{\text{out}}}`) dielectric constants combined in :math:`{\beta=\frac{\epsilon_{\text{in}}}{\epsilon_{\text{out}}}}`, atomic partial charges :math:`{q_{A/B}}`, and the electrostatic size of the solute :math:`{\mathcal{A}_{\text{det}}}`. \ :footcite:`ehlert2021`
:math:`{f(R_{AB}, a_{A}, a_{B})}` is the interaction kernel with the Born radii :math:`{a_{A/B}}` and can take two forms, either 

.. math::
   f(R_{AB}, a_{A}, a_{B}) = \left( R_{AB}^2 + a_{A} a_{B} \exp\left[-\frac{R_{AB}^2}{4 a_{A} a_{B}} \right] \right)^{\frac{1}{2}}

proposed by Still (default in GBSA), or the more recent P16 kernel (default for ALPB): 

.. math::
   f(R_{AB}, a_{A}, a_{B}) = R_{AB} + \sqrt{a_{A} a_{B}} \left(1+\frac{1.028 R_{AB}}{16 \sqrt{a_{A} a_{B}}} \right)^{16}

For specific polar interactions, an atom-wise hydrogen bonding correction is introduced:

.. math::
   \Delta G^{\text{HB}}_{\text{polar}} = \sum \Delta G^{\text{HB}}_{\text{A}}

In addition to the polar contribution, the non-polar contribution is included with a cavity dispersion solvation term (CDS) based on the atomic surface tension :math:`\gamma_{A}` and the solvent-accessible surface area (SASA) :math:`\sigma_{A}`: 

.. math::
   \Delta G^{\text{CDS}}_{\text{npol}} = \sum_{A} \gamma_{A} \sigma_{A}

An additional empirical constant shift is applied to the solvation free energy.
A solution state correction can be activated but is not included by default.

GBSA solvation model
--------------------

The generalized Born solvation model (GBSA) is a simplified version of ALPB in the limit of an ideal conductor environment (:math:`{\epsilon_{\text{out}}}\rightarrow \infty` and :math:`{\beta\rightarrow 0}`).
As for ALPB, CDS and a constant shift shift are applied, while a solution state correction can be activated (only if the solvent is specified by name).

Domain decomposition continuum solvation models
-----------------------------------------------

The domain decomposition solvation models implemented through the ddX library provide efficient formulations of the classical COSMO, CPCM, and PCM continuum models. The physical picture remains unchanged here, with the solute being embedded in a molecular cavity :math:`\Omega` that is surrounded by a homogeneous polarizable continuum. The solvent response is represented by an apparent surface charge distribution :math:`\sigma` on the cavity boundary :math:`\Gamma`, which generates the reaction field potential :math:`\phi^\sigma`. 

In the COSMO model, the solvent is first approximated by an ideal conductor (:math:`\epsilon \rightarrow \infty`), such that the total electrostatic potential must vanish at the cavity boundary. Therefore, the reaction field potential exactly cancels the molecular electrostatic potential :math:`\phi^{\rho_0}` on :math:`\Gamma`. Inside the cavity, the reaction field potential is harmonic, such that the following boundary value condition is fulfilled.

.. math::

   -\nabla^2 \phi^\sigma(\mathbf{r}) = 0, \qquad \mathbf{r} \in \Omega ,
   \phi^\sigma(\mathbf{s}) = -\phi^{\rho_0}(\mathbf{s}), \qquad \mathbf{s} \in \Gamma

The central idea of the domain decomposition formulation is to rewrite this global boundary value problem as a set of coupled local problems. For this purpose, the molecular cavity is decomposed back into its initial union of atom-centered spheres,

.. math::

   \Omega = \bigcup_{A=1}^{N_\mathrm{at}} \Omega_A .

The surface of each atomic sphere is further separated into an exposed (e) and buried (interior, i) part,

:: math::

   \Gamma_A = \Gamma_A^\mathrm{e} \cup \Gamma_A^\mathrm{i} .

Thus, the boundary value problem can be rewritten as a set of :math:`N_\mathrm{at}` coupled sphere-wise boundary value problems,

.. math::

   -\nabla^2\phi_A(\mathbf{s}) = 0, \qquad \mathbf{r} \in \Omega_A ,
   \phi_A(\mathbf{s}) = -\phi^{\rho_0}(\mathbf{s}), \qquad \mathbf{s} \in \Gamma_A^\mathrm{e} ,
   \phi_A(\mathbf{s}) = \frac{1}{|\mathcal{N}(A,\mathbf{s})|} \sum_{B \in \mathcal{N}(A,\mathbf{s})} \phi_B(\mathbf{s}), \qquad \mathbf{s} \in \Gamma_A^\mathrm{i},

where the potential on the buried surface is obtained from the potentials of all intersecting spheres at the considered surface point. Here, :math:`\mathcal{N}(A,\mathbf{s})` defines the list of neighbouring speres of sphere :math:`A` at surface point :math:`\mathbf{s}`. These local problems can be expressed in the known integral form by representing the reaction field potential on each sphere through a single-layer potential,

.. math::

   \phi_A(\mathbf{s}) = \int_{\Gamma_A} \frac{\sigma_A(\mathbf{s}')}{|\mathbf{s}-\mathbf{s}'|}
   \,d\mathbf{s}' .

The local surface charge distributions :math:`\sigma_A(\mathrm{s})` are then expanded in real spherical harmonics,

.. math::

   \sigma_A(\mathbf{s}) = \sum_{\ell=0}^{\ell_\mathrm{max}} \sum_{m=-\ell}^{\ell} \[X_A\]_{\ell m} Y_{\ell m}(\mathbf{s}) ,

and integration is performed numerically on a unit sphere via Lebedev quadratures. This finally results in a linear system

.. math::

   \mathbf{L} \mathbf{X} = \mathbf{g} ,

where the vector :math:`\mathbf{X}` contains the expansion coefficients of the local apparent surface charges, :math:`\mathbf{g}` is the corresponding spherical harmonics expansion of the free-space potential generated by the solutes's charge distribution, and :math:`\mathbf{L}` is a naturally block-sparse coupling matrix.
This system is solved iteratively in a very efficient manner, for instance by Jacobi or DIIS procedures, which allows to compute the solvation contribution with linear-scaling behavior.

In xTB, the solute charge density is represented by atom-centered point charges. The final solvation free energy is therefore obtained through a dot product of the vector containing the partial charges, :math:`\mathbf{\Psi}`, and the expansion coefficients,

.. math::

   \Delta G_\mathrm{polar} = \frac{1}{2}f(\epsilon,x} \langle \mathbf{X}, \mathbf{\Psi} \rangle ,

where the factor

.. math::

   f(\epsilon,x) = \frac{\epsilon -1}{\epsilon +x} 

scales the surface charges, or equivalently the energy, to account for the fact that the real dielectric permittivity of the solvent is in fact finite. As established, it further distinguishes between COSMO (:math:`x=0.5`) and CPCM (:math:`x=0`).

In PCM, the continuum is not treated as an ideal conductor, leading to both single- and double-layer potential contributions. Adressing the general integral equation formalism (IEF) of PCM by the same sphere-wise framework used for ddCOSMO/ddCPCM, the ddPCM equation can be written as 

.. math::

   \hat{\mathcal{R}}_\epsilon \hat{\mathcal{S}}\sigma = -\hat{\mathcal{R}}_\infty \phi^{\rho_0} ,

with the operators

.. math::

   \hat{\mathcal{R}}_\epsilon = \frac{2\pi}{f(\epsilon)}\hat{1} -\hat{\mathcal{D}}

   \hat{\mathcal{R}}_\infty = 2\pi\hat{1}-\hat{\mathcal{D}} ,

and the single- and double-layer operators :math:`\hat{\mathcal{S}` and `\hat{\mathcal{D}`. The ddPCM model solves this problem by introducing an intermediate potential :math:`\phi_\epsilon`, which splits the problem into 

.. math::

   \hat{\mathcal{R}}_\epsilon \phi_\epsilon = \hat{\mathcal{R}}_\infty \phi^{\rho_0}

and 

.. math::

   \hat{\mathcal{S}}\sigma = -\phi_\epsilon .

While the second equation has the same structure as the single-layer problem solved in ddCOSMO/ddCPCM, the first equation stands out as the additional PCM problem.

For further reading, the user is referred to comprehensive reviews \ :footcite:`stamm_how_2019` \ :footcite:`nottoli_ddx_2024` and the references cited therein.


CPCM solvation model 
--------------------

The conductor-like polarizable continuum solvation model is implemented based on the domain-decomposition approach and is currently available only for the polar part :math:`{\Delta G_{\text{polar}}}`.

Solution state correction
-------------------------

For solvation free energies, the state of the inital gas and final liquid solution can be changed with a solution state correction.
By default no solution state correction is applied (gsolv, default), which is comparable with most other solvation models (SMD, COSMO-RS, ...).
For normal production runs, the option ``bar1mol`` should be used. For explicit comparisons with ``reference`` state corrected COSMO-RS, the ``reference`` option should be used (includes solvent-specific correction for infinite dilution).
Solution state correction is available for the ALPB and GBSA solvation models.

================== ====================================================================
 Name               Definition
================== ====================================================================
 gsolv (default)    1 L of ideal gas and 1 L of liquid solution 
 bar1mol            1 bar of ideal gas and 1 mol/L liquid solution 
 reference          1 bar of ideal gas and 1 mol/L liquid solution at infinite dilution
================== ====================================================================


Literature
----------

.. footbibliography::
