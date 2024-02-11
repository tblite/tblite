"""
"""

from typing import Any, Callable, Generator, Optional, Sequence, Tuple

import numpy as np

try:
    import berny
except ModuleNotFoundError as e:
    raise ModuleNotFoundError("This submodule requires pyberny installed") from e

from .interface import Calculator, Result
from .utils import to_number

BernyYield = Optional[Tuple[float, np.ndarray]]
BernySend = Tuple[Sequence[tuple[str, Sequence[float]]], Optional[Sequence[float]]]
Callback = Callable[
    [np.ndarray, np.ndarray, Optional[np.ndarray], float, np.ndarray, Result], None
]


def Solver(
    method: str,
    charge: Optional[float] = None,
    uhf: Optional[int] = None,
    *,
    callback: Callback = None,
    parameters: dict[str, Any] = {},
    context_options: dict[str, Any] = {},
) -> Generator[BernyYield, BernySend, None]:
    atoms, lattice = yield
    numbers = np.asarray([to_number(sp) for sp, _ in atoms])
    positions = np.asarray([coord for _, coord in atoms]) * berny.coords.angstrom
    if lattice is not None:
        lattice = np.asarray(lattice) * berny.coords.angstrom
    calc = Calculator(
        method=method,
        numbers=numbers,
        positions=positions,
        charge=charge,
        uhf=uhf,
        lattice=lattice,
        **context_options,
    )

    options, interactions = calc.check_parameters(parameters)
    for option, value in options.items():
        calc.set(option, value)
    for interaction, value in interactions.items():
        calc.add(interaction, *value)

    if callback is not None:
        callback(numbers, positions, lattice, None, None, None)
    res = None
    while True:
        positions = np.asarray([coord for _, coord in atoms]) * berny.coords.angstrom
        if lattice is not None:
            lattice = np.asarray(lattice) * berny.coords.angstrom
        calc.update(positions=positions, lattice=lattice)
        res = calc.singlepoint(res)
        energy, gradient = res.get("energy"), res.get("gradient")
        if callback is not None:
            callback(numbers, positions, lattice, energy, gradient, res)
        atoms, lattice = yield energy, gradient / berny.coords.angstrom
