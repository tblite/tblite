# This file is part of tblite.
# SPDX-Identifier: LGPL-3.0-or-later
#
# tblite is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tblite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with tblite.  If not, see <https://www.gnu.org/licenses/>.
"""
The Python API of *tblite* natively support integration with the atomic simulation environment (`ASE`_).
By constructing a calculator most functionality of ASE is readily available.
For details on building the Python API checkout the :ref:`installation guide <python-build>`.

.. _ase: https://wiki.fysik.dtu.dk/ase/
"""

try:
    import ase
except ModuleNotFoundError:
    raise ModuleNotFoundError("This submodule requires ASE installed")


from typing import List, Optional

from .interface import Calculator
import ase.calculators.calculator
from ase.atoms import Atoms
from ase.units import Hartree, Bohr, kB


class TBLite(ase.calculators.calculator.Calculator):
    """
    ASE calculator for using xTB Hamiltonians from the tblite library.
    Supported properties by this calculator are:

    - energy (free_energy)
    - forces
    - stress
    - dipole
    - charges

    Supported keywords are

    ======================== ================= ============================================
     Keyword                  Default           Description
    ======================== ================= ============================================
     method                   "GFN2-xTB"        Underlying method for energy and forces
     accuracy                 1.0               Numerical accuracy of the calculation
     electronic_temperature   300.0             Electronic temperatur in Kelvin
     max_iterations           250               Iterations for self-consistent evaluation
     cache_api                True              Reuse generate API objects (recommended)
    ======================== ================= ============================================

    Example
    -------

    An ASE calculator can be constructed by using the *TBLite* class provided by the *tblite.ase* module.
    For example to perform a single point calculation for a CO\ :sub:`2` crystal use

    >>> from tblite.ase import TBLite
    >>> from ase.atoms import Atoms
    >>> import numpy as np
    >>> atoms = Atoms(
    ...     symbols="C4O8",
    ...     positions=np.array(
    ...         [
    ...             [0.9441259872, 0.9437851680, 0.9543505632],
    ...             [3.7179966528, 0.9556570368, 3.7316862240],
    ...             [3.7159517376, 3.7149292800, 0.9692330016],
    ...             [0.9529872864, 3.7220864832, 3.7296981120],
    ...             [1.6213905408, 1.6190616096, 1.6313879040],
    ...             [0.2656685664, 0.2694175776, 0.2776540416],
    ...             [4.3914553920, 1.6346256864, 3.0545920000],
    ...             [3.0440834880, 0.2764611744, 4.4080419264],
    ...             [4.3910577696, 3.0416409504, 0.2881058304],
    ...             [3.0399936576, 4.3879335936, 1.6497353376],
    ...             [0.2741322432, 4.4003734944, 3.0573754368],
    ...             [1.6312174944, 3.0434586528, 4.4023048032],
    ...         ]
    ...     ),
    ...     cell=np.array([5.68032, 5.68032, 5.68032]),
    ...     pbc=np.array([True, True, True]),
    ... )
    >>> atoms.calc = TBLite(method="GFN1-xTB")
    >>> atoms.get_potential_energy()  # result in eV
    -1257.0943962462964

    The resulting calculator can be used like most ASE calculator, *e.g.* for optimizing geometries.
    """

    implemented_properties = [
        "energy",
        "forces",
        "charges",
        "dipole",
        "stress",
    ]

    default_parameters = {
        "method": "GFN2-xTB",
        "accuracy": 1.0,
        "max_iterations": 250,
        "electronic_temperature": 300.0,
        "cache_api": True,
    }

    _res = None
    _xtb = None

    def __init__(
        self,
        atoms: Optional[Atoms] = None,
        **kwargs,
    ):
        """
        Construct the TBLite base calculator object.
        """

        ase.calculators.calculator.Calculator.__init__(self, atoms=atoms, **kwargs)

    def set(self, **kwargs) -> dict:
        """
        Set new parameters to TBLite. Will automatically reconstruct the underlying
        model in case critical parameters change.

        Example
        -------
        >>> from ase.build import molecule
        >>> from tblite.ase import TBLite
        >>> atoms = molecule("H2O")
        >>> atoms.calc = TBLite(method="GFN2-xTB")
        >>> atoms.get_potential_energy()
        -137.96777625229421
        >>> atoms.calc.set(method="GFN1-xTB")
        {'method': 'GFN1-xTB'}
        >>> atoms.get_potential_energy()
        -156.9675057724589
        """

        changed_parameters = ase.calculators.calculator.Calculator.set(self, **kwargs)

        # Always reset the calculation if parameters change
        if changed_parameters:
            self.reset()

        # If the method is changed, invalidate the cached calculator as well
        if "method" in changed_parameters:
            self._xtb = None
            self._res = None

        # Minor changes can be updated in the API calculator directly
        if self._xtb is not None:
            if "accuracy" in changed_parameters:
                self._xtb.set("accuracy", self.parameters.accuracy)

            if "electronic_temperature" in changed_parameters:
                self._xtb.set(
                    "temperature", self.parameters.electronic_temperature * kB / Hartree
                )

            if "max_iterations" in changed_parameters:
                self._xtb.set("max-iter", self.parameters.max_iterations)

        return changed_parameters

    def reset(self) -> None:
        """
        Clear all information from old calculation. This will only remove the cached
        API objects in case the `cache_api` is set to False.
        """
        ase.calculators.calculator.Calculator.reset(self)

        if not self.parameters.cache_api:
            self._xtb = None
            self._res = None

    def _check_api_calculator(self, system_changes: List[str]) -> None:
        """Check state of API calculator and reset if necessary"""

        # Changes in positions and cell parameters can use a normal update
        _reset = system_changes.copy()
        if "positions" in _reset:
            _reset.remove("positions")
        if "cell" in _reset:
            _reset.remove("cell")

        # Invalidate cached calculator and results object
        if _reset:
            self._xtb = None
            self._res = None
        else:
            if system_changes and self._xtb is not None:
                try:
                    _cell = self.atoms.cell
                    self._xtb.update(
                        self.atoms.positions / Bohr,
                        _cell / Bohr,
                    )
                # An exception in this part means the geometry is bad,
                # still we will give a complete reset a try as well
                except RuntimeError:
                    self._xtb = None
                    self._res = None

    def _create_api_calculator(self) -> Calculator:
        """Create a new API calculator object"""

        try:
            _cell = self.atoms.cell
            _periodic = self.atoms.pbc
            _charge = self.atoms.get_initial_charges().sum()
            _uhf = int(self.atoms.get_initial_magnetic_moments().sum().round())

            calc = Calculator(
                self.parameters.method,
                self.atoms.numbers,
                self.atoms.positions / Bohr,
                _charge,
                _uhf,
                _cell / Bohr,
                _periodic,
            )
            calc.set("accuracy", self.parameters.accuracy)
            calc.set(
                "temperature", self.parameters.electronic_temperature * kB / Hartree
            )
            calc.set("max-iter", self.parameters.max_iterations)

        except RuntimeError:
            raise ase.calculators.calculator.InputError(
                "Cannot construct calculator for TBLite"
            )

        return calc

    def calculate(
        self,
        atoms: Optional[Atoms] = None,
        properties: List[str] = None,
        system_changes: List[str] = ase.calculators.calculator.all_changes,
    ) -> None:
        """
        Perform actual calculation with by calling the TBLite API

        Example
        -------

        >>> from ase.build import molecule
        >>> from tblite.ase import TBLite
        >>> calc = TBLite(method="GFN2-xTB")
        >>> calc.calculate(molecule("H2O"))
        >>> calc.get_potential_energy()
        -137.96777625229421
        >>> calc.calculate(molecule("CH4"))
        >>> calc.get_potential_energy()
        -113.60956621093894

        Raises
        ------
        ase.calculators.calculator.InputError
            on invalid input passed to the interface module

        ase.calculators.calculator.CalculationFailed
            in case of an `RuntimeError` in the library
        """

        if not properties:
            properties = ["energy"]
        ase.calculators.calculator.Calculator.calculate(
            self, atoms, properties, system_changes
        )

        self._check_api_calculator(system_changes)

        if self._xtb is None:
            self._xtb = self._create_api_calculator()

        try:
            self._res = self._xtb.singlepoint(self._res)
        except RuntimeError:
            raise ase.calculators.calculator.CalculationFailed(
                "TBLite could not evaluate input"
            )

        # These properties are garanteed to exist for all implemented calculators
        self.results["energy"] = self._res.get("energy") * Hartree
        self.results["free_energy"] = self.results["energy"]
        self.results["forces"] = -self._res.get("gradient") * Hartree / Bohr
        self.results["charges"] = self._res.get("charges")
        self.results["dipole"] = self._res.get("dipole") * Bohr
        # stress tensor is only returned for periodic systems
        if self.atoms.pbc.any():
            _stress = self._res.get("virial") * Hartree / self.atoms.get_volume()
            self.results["stress"] = _stress.flat[[0, 4, 8, 5, 2, 1]]
