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
`ASE calculator <https://wiki.fysik.dtu.dk/ase/>`_ implementation for the tblite library.

Supported properties by this calculator are:

- energy (free_energy)
- forces
- stress
- dipole
- charges

Supported keywords are

======================== ============ ============================================
 Keyword                  Default      Description
======================== ============ ============================================
 method                   "GFN2-xTB"   Underlying method for energy and forces
 accuracy                 1.0          Numerical accuracy of the calculation
 electronic_temperature   300.0        Electronic temperatur in Kelvin
 max_iterations           250          Iterations for self-consistent evaluation
 cache_api                True         Reuse generate API objects (recommended)
======================== ============ ============================================
"""

try:
    import ase
except ModuleNotFoundError:
    raise ModuleNotFoundError("This submodule requires ASE installed")


from typing import List, Optional

from .interface import Calculator, Result
import ase.calculators.calculator
from ase.atoms import Atoms
from ase.units import Hartree, Bohr, kB


class TBLite(ase.calculators.calculator.Calculator):
    """
    ASE calculator for using xTB Hamiltonians from the tblite library.
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
        """Construct the TBLite base calculator object."""

        ase.calculators.calculator.Calculator.__init__(self, atoms=atoms, **kwargs)

    def set(self, **kwargs) -> dict:
        """Set new parameters to TBLite"""

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
        """Clear all information from old calculation"""
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
        """Perform actual calculation with by calling the TBLite API"""

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
