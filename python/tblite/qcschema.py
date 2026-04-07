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
Integration with the `QCArchive infrastructure <https://qcarchive.molssi.org>`_.

This module provides a function to run QCSchema input or QCElemental Atomic Input
through the TBLite calculator and return the result as a QCElemental Atomic Result.

If the QCElemental package is installed the ``tblite.qcschema`` module becomes
importable and provides the ``run_schema`` function supporting QCSchema v1.
If the QCElemental package is >=0.50.0, ``tblite.qcschema`` supports QCSchema v1
and v2, returning whichever version was submitted. Note that Python 3.14+ only
works with QCSchema v2 due to Pydantic restrictions.

The model support the following methods:

- **GFN2-xTB**: 
  Self-consistent extended tight binding Hamiltonian with
  anisotropic second order electrostatic contributions,
  third order on-site contributions and self-consistent D4 dispersion.
  Geometry, frequency and non-covalent interactions parametrisation for
  elements up to Z=86.

- **GFN1-xTB**:
  Self-consistent extended tight binding Hamiltonian with
  isotropic second order electrostatic contributions and
  third order on-site contributions.
  Geometry, frequency and non-covalent interactions parametrisation for
  elements up to Z=86.

- **IPEA1-xTB**:
  Special parametrisation for the GFN1-xTB Hamiltonian to improve the
  description of vertical ionisation potentials and electron affinities.
  Uses additional diffuse s-functions on light main group elements.
  Parametrised up to Z=86.

Supported keywords are:

=================== ==================================== =========================================
 name                description                          default
=================== ==================================== =========================================
 accuracy            Numerical thresholds for SCC         float (1.0)
 guess               Initial guess for wavefunction       integer (0 == SAD)
 max-iter            Maximum number of SCC iterations     integer (250)
 mixer-damping       Parameter for the SCC mixer          float (0.4)
 save-integrals      Keep integral matrices in results    0 (False)
 temperature         Electronic temperature for filling   float (9.500e-4)
 verbosity           Set verbosity of printout            integer (1)
 electric-field      Uniform electric field               Field vector
 spin-polarization   Spin polarization                    Scaling factor
 alpb-solvation      ALPB implicit solvation              Solvent name, solution state (optional)
 gbsa-solvation      GBSA implicit solvation              Solvent name, solution state (optional)
 cpcm-solvation      CPCM implicit solvation              Epsilon
 gbe-solvation       GBε implicit solvation               Epsilon, Born kernel
 gb-solvation        GB implicit solvation                Epsilon, Born kernel
=================== ==================================== =========================================
"""

import sys
from io import StringIO
from typing import Any, Dict, Literal, overload, Union

import numpy as np

from .exceptions import TBLiteRuntimeError, TBLiteTypeError, TBLiteValueError
from .interface import Calculator
from .library import get_version

if sys.version_info < (3, 14):
    try:
        import qcelemental.models.v1 as qcel_v1
    except ModuleNotFoundError:
        import qcelemental.models as qcel_v1
else:
    qcel_v1 = None

try:
    import qcelemental.models.v2 as qcel_v2
except ModuleNotFoundError:
    qcel_v2 = None


if qcel_v1 is None and qcel_v2 is None:
    raise ModuleNotFoundError(
        "The qcelemental package is required for qcschema support. "
        "Please install it with 'pip install qcelemental'."
    )

SUPPORTED_DRIVERS = {
    "energy",
    "gradient",
    # "hessian",
    "properties",
}


if qcel_v1 is not None:

    @overload
    def get_provenance(schema_version: Literal[1]) -> "qcel_v1.Provenance": ...

    @overload
    def get_error(
        input_data: "qcel_v1.AtomicInput",
        error: Union[Dict[str, Any], "qcel_v1.ComputeError"],
        schema_version: int,
    ) -> "qcel_v1.AtomicResult": ...

    @overload
    def run_schema(
        input_data: Union[Dict[str, Any], "qcel_v1.AtomicInput"],
    ) -> Union["qcel_v1.AtomicResult", "qcel_v1.FailedOperation"]: ...


if qcel_v2 is not None:

    @overload
    def get_provenance(schema_version: Literal[2]) -> "qcel_v2.Provenance": ...

    @overload
    def get_error(
        input_data: "qcel_v2.AtomicInput",
        error: Union[Dict[str, Any], "qcel_v2.ComputeError"],
        schema_version: int,
    ) -> "qcel_v2.FailedOperation": ...

    @overload
    def run_schema(
        input_data: Union[Dict[str, Any], "qcel_v2.AtomicInput"],
    ) -> Union["qcel_v2.AtomicResult", "qcel_v2.FailedOperation"]: ...


def get_provenance(schema_version: Literal[1, 2]):
    """
    Returns a QCSchema provenance model.

    Returns
    -------
    qcel.models.Provenance
        A QCSchema provenance model.
    """
    prov = dict(
        creator="tblite",
        version=".".join([str(v) for v in get_version()]),
        routine="tblite.qcschema.run_schema",
    )
    if schema_version == 1:
        return qcel_v1.Provenance(**prov)
    elif schema_version == 2:
        return qcel_v2.Provenance(**prov)
    else:
        raise ValueError(
            f"Unsupported QCSchema version: {schema_version}. Only v1 and v2 are supported."
        )


def get_error(input_data, error, schema_version):
    if schema_version == 1:
        if not isinstance(error, qcel_v1.ComputeError):
            error = qcel_v1.ComputeError(**error)

        if not isinstance(input_data, qcel_v1.AtomicInput):
            input_data = qcel_v1.AtomicInput(**input_data)

        return_data = input_data.dict()
        return_data.update(
            error=error,
            success=False,
            return_result={
                "energy": 0.0,
                "gradient": np.zeros(input_data.molecule.geometry.shape),
                "hessian": np.zeros(
                    (
                        input_data.molecule.geometry.size,
                        input_data.molecule.geometry.size,
                    )
                ),
                "properties": {},
            }[input_data.driver],
            properties={},
            provenance=get_provenance(schema_version),
        )
        return qcel_v1.AtomicResult(**return_data)

    elif schema_version == 2:
        if not isinstance(error, qcel_v2.ComputeError):
            error = qcel_v2.ComputeError(**error)

        return qcel_v2.FailedOperation(input_data=input_data, error=error)

    else:
        raise ValueError(
            f"Unsupported QCSchema version: {schema_version}. Only v1 and v2 are supported."
        )


class _Logger:
    def __init__(self):
        self._buffer = StringIO()

    def __call__(self, message: str) -> None:
        print(message)
        self._buffer.write(message + "\n")

    def __str__(self) -> str:
        return self._buffer.getvalue()


def run_schema(input_data):
    """Runs a QCSchema input through the QCEngine stack and returns a QCSchema result.

    Parameters
    ----------
    input_data : qcel.models.AtomicInput
        A QCSchema v1 or v2 input dictionary or model.

    Returns
    -------
    qcel.models.AtomicResult
        A QCSchema v1 or v2 result model.
    """
    if qcel_v2 is not None and isinstance(input_data, qcel_v2.AtomicInput):
        atomic_input = input_data
    elif qcel_v1 is not None and isinstance(input_data, qcel_v1.AtomicInput):
        atomic_input = input_data
    elif qcel_v2 is not None and input_data.get("specification"):
        atomic_input = qcel_v2.AtomicInput(**input_data)
    elif qcel_v1 is not None:
        atomic_input = qcel_v1.AtomicInput(**input_data)
    else:
        raise ValueError(
            "Input data is not a valid QCSchema AtomicInput for either v1 or v2."
        )

    schema_version = atomic_input.schema_version
    if schema_version == 1:
        ret_data = atomic_input.dict()
        input_keywords = atomic_input.keywords
        input_method = atomic_input.model.method
        input_basis = atomic_input.model.basis
        input_driver = atomic_input.driver
    elif schema_version == 2:
        ret_data = {
            "input_data": atomic_input,
            "extras": {},
            "molecule": atomic_input.molecule,
        }
        input_keywords = atomic_input.specification.keywords
        input_method = atomic_input.specification.model.method
        input_basis = atomic_input.specification.model.basis
        input_driver = atomic_input.specification.driver
    else:
        raise ValueError(
            f"Unsupported QCSchema version: {schema_version}. Only v1 and v2 are supported."
        )

    if input_driver not in SUPPORTED_DRIVERS:
        driver_name = (
            input_driver.name
            if hasattr(input_driver, "name")
            else str(input_driver)
        )
        return get_error(
            input_data,
            dict(
                error_type="input_error",
                error_message=f"Driver '{driver_name}' is not supported by tblite.",
            ),
            schema_version,
        )

    if input_method not in Calculator._loader:
        return get_error(
            input_data,
            dict(
                error_type="input_error",
                error_message=f"Model '{input_method}' is not supported by tblite.",
            ),
            schema_version,
        )

    if input_basis is not None:
        return get_error(
            input_data,
            dict(
                error_type="input_error",
                error_message="Basis sets are not supported by tblite.",
            ),
            schema_version,
        )

    keywords = {
        key: value
        for key, value in input_keywords.items()
        if key in Calculator._setter
    }
    interaction = {
        key: value
        for key, value in input_keywords.items()
        if key in Calculator._interaction
    }
    unknown_keywords = (
        set(input_keywords) - set(keywords) - set(interaction)
    )
    if unknown_keywords:
        return get_error(
            input_data,
            dict(
                error_type="input_error",
                error_message=f"Unknown keywords: {', '.join(unknown_keywords)}.",
            ),
            schema_version,
        )

    logger = _Logger()
    try:
        calc = Calculator(
            method=input_method,
            numbers=atomic_input.molecule.atomic_numbers,
            positions=atomic_input.molecule.geometry,
            charge=atomic_input.molecule.molecular_charge,
            uhf=atomic_input.molecule.molecular_multiplicity - 1,
            color=False,
            logger=logger,
        )
        for key, value in keywords.items():
            calc.set(key, value)

        for key, value in interaction.items():
            if not isinstance(value, list):
                value = [value]
            calc.add(key, *value)

        result = calc.singlepoint().dict()

        if schema_version == 1:
            AtProp = qcel_v1.AtomicResultProperties
        elif schema_version == 2:
            AtProp = qcel_v2.AtomicProperties

        properties = AtProp(
            return_energy=result["energy"],
            return_gradient=result["gradient"],
            calcinfo_natom=result["natoms"],
            calcinfo_nbasis=result["norbitals"],
            calcinfo_nmo=result["norbitals"],
            scf_dipole_moment=result["dipole"],
            scf_quadrupole_moment=result["quadrupole"][
                [0, 1, 3, 1, 2, 4, 3, 4, 5]
            ],
            scf_total_energy=result["energy"],
            scf_total_gradient=result["gradient"],
        )

        return_result = {
            "energy": result["energy"],
            "gradient": result["gradient"],
            "properties": {
                "dipole": result["dipole"],
                "mulliken_charges": result["charges"],
                "mayer_indices": result["bond-orders"],
            },
        }[input_driver]

    except (TBLiteTypeError, TBLiteValueError) as e:
        return get_error(
            input_data,
            error = dict(
                error_type="input_error",
                error_message=str(e),
            ),
            schema_version=schema_version,
        )

    except TBLiteRuntimeError as e:
        return get_error(
            input_data,
            error = dict(
                error_type="execution_error",
                error_message=str(e),
            ),
            schema_version=schema_version,
        )

    ret_data.update(
        success=True,
        return_result=return_result,
        properties=properties,
        provenance=get_provenance(schema_version),
        stdout=str(logger),
    )

    if schema_version == 1:
        return qcel_v1.AtomicResult(**ret_data)

    return qcel_v2.AtomicResult(**ret_data)
