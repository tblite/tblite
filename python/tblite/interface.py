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
Definition of the basic interface to library for most further integration in
other Python frameworks. The classes defined here allow a more Pythonic usage
of the library in actual workflows than the low-level access provided in the
CFFI generated wrappers.
"""

import numpy as np
from typing import Optional

from . import library


class Structure:
    """
    .. Molecular structure data

    Represents a wrapped structure object in ``tblite``.
    The molecular structure data object has a fixed number of atoms
    and immutable atomic identifiers.
    """

    _mol = library.ffi.NULL

    def __init__(
        self,
        numbers: np.ndarray,
        positions: np.ndarray,
        charge: Optional[float] = None,
        uhf: Optional[int] = None,
        lattice: Optional[np.ndarray] = None,
        periodic: Optional[np.ndarray] = None,
    ):
        """Create new molecular structure data"""
        if positions.size % 3 != 0:
            raise ValueError("Expected tripels of cartesian coordinates")

        if 3 * numbers.size != positions.size:
            raise ValueError("Dimension missmatch between numbers and positions")

        self._natoms = len(numbers)
        _numbers = np.ascontiguousarray(numbers, dtype="i4")
        _positions = np.ascontiguousarray(positions, dtype=float)

        _charge = _ref("double", charge)
        _uhf = _ref("int", uhf)

        if lattice is not None:
            if lattice.size != 9:
                raise ValueError("Invalid lattice provided")
            _lattice = np.ascontiguousarray(lattice, dtype="float")
        else:
            _lattice = None

        if periodic is not None:
            if periodic.size != 3:
                raise ValueError("Invalid periodicity provided")
            _periodic = np.ascontiguousarray(periodic, dtype="bool")
        else:
            _periodic = None

        self._mol = library.new_structure(
            self._natoms,
            _cast("int*", _numbers),
            _cast("double*", _positions),
            _charge,
            _uhf,
            _cast("double*", _lattice),
            _cast("bool*", _periodic),
        )

    def __len__(self):
        return self._natoms

    def update(
        self,
        positions: np.ndarray,
        lattice: Optional[np.ndarray] = None,
    ) -> None:
        """Update coordinates and lattice parameters, both provided in
        atomic units (Bohr).
        The lattice update is optional also for periodic structures.

        Generally, only the cartesian coordinates and the lattice parameters
        can be updated, every other modification, regarding total charge,
        total spin, boundary condition, atomic types or number of atoms
        requires the complete reconstruction of the object.
        """

        if 3 * len(self) != positions.size:
            raise ValueError("Dimension missmatch for positions")
        _positions = np.ascontiguousarray(positions, dtype="float")

        if lattice is not None:
            if lattice.size != 9:
                raise ValueError("Invalid lattice provided")
            _lattice = np.ascontiguousarray(lattice, dtype="float")
        else:
            _lattice = None

        library.update_structure_geometry(
            self._mol,
            _cast("double*", _positions),
            _cast("double*", _lattice),
        )


class Result:
    """
    .. Calculation result and restart data
    """

    _res = library.ffi.NULL
    _natoms = 0
    _getter = {
        "energy": library.get_energy,
        "gradient": library.get_gradient,
        "virial": library.get_virial,
        "charges": library.get_charges,
        "dipole": library.get_dipole,
        "quadrupole": library.get_quadrupole,
        "orbital-energies": library.get_orbital_energies,
        "orbital-occupations": library.get_orbital_occupations,
        "orbital-coefficients": library.get_orbital_coefficients,
    }
    _setter = {}

    def __init__(self, other=None):
        if other is not None:
            self._res = library.copy_result(other._res)
        else:
            self._res = library.new_result()

    def get(self, attribute: str):
        """Get a quantity stored instade the result container"""

        if attribute not in self._getter:
            raise ValueError(f"Attribute '{attribute}' is not available in this result")

        return self._getter[attribute](self._res)

    def set(self, attribute: str, value):
        """Get a quantity stored instade the result container"""

        if attribute not in self._setter:
            raise ValueError(f"Attribute '{attribute}' cannot be set in this result")

        self._setter[attribute](self._res, value)

    def dict(self) -> dict:
        """Return all quantities inside the result container as dict"""
        res = {}

        for key in self._getter:
            try:
                res[key] = self.get(key)
            except RuntimeError:
                pass

        return res


class Calculator(Structure):
    """
    .. Singlepoint calculator
    """

    _ctx = library.ffi.NULL
    _calc = library.ffi.NULL
    _loader = {
        "GFN2-xTB": library.new_gfn2_calculator,
        "GFN1-xTB": library.new_gfn1_calculator,
        "IPEA1-xTB": library.new_ipea1_calculator,
    }
    _setter = {
        "max-iter": library.set_calculator_max_iter,
        "accuracy": library.set_calculator_accuracy,
        "mixer-damping": library.set_calculator_mixer_damping,
        "temperature": library.set_calculator_temperature,
    }

    def __init__(
        self,
        method: str,
        numbers: np.ndarray,
        positions: np.ndarray,
        charge: Optional[float] = None,
        uhf: Optional[int] = None,
        lattice: Optional[np.ndarray] = None,
        periodic: Optional[np.ndarray] = None,
    ):
        """Construct new calculator object for a given structure"""
        Structure.__init__(self, numbers, positions, charge, uhf, lattice, periodic)

        self._ctx = library.new_context()
        if method not in self._loader:
            raise ValueError(f"Method '{method}' is not available for this calculator")
        self._calc = self._loader[method](self._ctx, self._mol)

    def set(self, attribute: str, value) -> None:
        """Set an attribute in the calculator instance"""

        if attribute not in self._setter:
            raise ValueError(
                f"Attribute '{attribute}' is not supported in this calculator"
            )
        self._setter[attribute](self._ctx, self._calc, value)

    def singlepoint(self, res: Optional[Result] = None, copy: bool = False) -> Result:
        """Perform actual single point calculation"""

        _res = Result(res) if copy or res is None else res

        library.get_singlepoint(self._ctx, self._mol, self._calc, _res._res)
        return _res


def _cast(ctype, array):
    """Cast a numpy array to an FFI pointer"""
    return (
        library.ffi.NULL
        if array is None
        else library.ffi.cast(ctype, array.ctypes.data)
    )


def _ref(ctype, value):
    """Create a reference to a value"""
    if value is None:
        return library.ffi.NULL
    ref = library.ffi.new(ctype + "*")
    ref[0] = value
    return ref
