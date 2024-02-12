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
Interface to the `pyberny`_ library for geometry optimization.

.. _pyberny: https://jan.hermann.name/pyberny/algorithm.html

Example
-------
>>> import numpy as np
>>> from berny import Berny, Geometry, optimize
>>> from tblite.berny import Solver
>>> geom = Geometry(
...     species=["O", "H", "H"],
...     coords=np.array([
...         [+0.000000, +0.000000, +0.119262],
...         [+0.000000, +0.763239, -0.477047],
...         [+0.000000, -0.763239, -0.477047],
...     ]),
... )
>>> opt = Berny(geom, maxsteps=10)
>>> solver = Solver("GFN2-xTB", parameters={"verbosity": 0})
>>> geom = optimize(opt, solver)
>>> opt.converged
True
>>> np.set_printoptions(suppress=True)
>>> geom.coords
array([[ 0.        , -0.        ,  0.10125553],
       [ 0.        ,  0.77171228, -0.46804376],
       [ 0.        , -0.77171228, -0.46804376]])
"""

from typing import Any, Callable, Dict, Generator, Optional, Sequence, Tuple

import numpy as np

try:
    import berny
except ModuleNotFoundError as e:
    raise ModuleNotFoundError("This submodule requires pyberny installed") from e

from .interface import Calculator, Result
from .utils import to_number

BernyYield = Optional[Tuple[float, np.ndarray]]
BernySend = Tuple[Sequence[Tuple[str, Sequence[float]]], Optional[Sequence[float]]]
Callback = Callable[
    [np.ndarray, np.ndarray, Optional[np.ndarray], float, np.ndarray, Result], None
]


def Solver(
    method: str,
    charge: Optional[float] = None,
    uhf: Optional[int] = None,
    *,
    callback: Callback = None,
    parameters: Dict[str, Any] = {},
    context_options: Dict[str, Any] = {},
) -> Generator[BernyYield, BernySend, None]:
    """
    Wraps a tblite calculator to make it a generator for use in optimizer.

    Parameters
    ----------
    method: str
        Name of the method to initialize the calculator with
    charge: float (optional)
        Total charge of the system
    uhf: float (optional)
        Number of unpaired electrons in the system
    callback: Callable (optional)
        Callback function used after each successful geometry step
    parameters: dict (optional)
        Options passed through to the singlepoint calculator
    context_options: dict (optional)
        Options passed through to the calculation context

    Example
    -------
    >>> import numpy as np
    >>> from berny import Geometry
    >>> from tblite.berny import Solver
    >>> geom = Geometry(
    ...     species=["O", "H", "H"],
    ...     coords=np.array([
    ...         [+0.000000, +0.000000, +0.119262],
    ...         [+0.000000, +0.763239, -0.477047],
    ...         [+0.000000, -0.763239, -0.477047],
    ...     ]),
    ... )
    >>> solver = Solver("GFN2-xTB", parameters={"verbosity": 0})
    >>> next(solver)
    >>> energy, gradient = solver.send((list(geom), geom.lattice))
    >>> energy
    array(-5.07022229)
    """
    atoms, lattice = yield
    numbers = np.asarray([to_number(sp) for sp, _ in atoms])
    positions = np.asarray([coord for _, coord in atoms]) * berny.angstrom
    if lattice is not None:
        lattice = np.asarray(lattice) * berny.angstrom
    calc = Calculator(
        method=method,
        numbers=numbers,
        positions=positions,
        charge=charge,
        uhf=uhf,
        lattice=lattice,
        **context_options,
    )

    options, interactions = calc.split_parameters(parameters)
    for option, value in options.items():
        calc.set(option, value)
    for interaction, value in interactions.items():
        calc.add(interaction, *value)

    res = None
    while True:
        positions = np.asarray([coord for _, coord in atoms]) * berny.angstrom
        if lattice is not None:
            lattice = np.asarray(lattice) * berny.angstrom
        calc.update(positions=positions, lattice=lattice)
        res = calc.singlepoint(res)
        energy, gradient = res.get("energy"), res.get("gradient")
        if callback is not None:
            callback(numbers, positions, lattice, energy, gradient, res)
        atoms, lattice = yield energy, gradient / berny.angstrom
