/* This file is part of tblite.
 * SPDX-Identifier: LGPL-3.0-or-later
 *
 * tblite is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tblite is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tblite.  If not, see <https://www.gnu.org/licenses/>.
**/

#pragma once

#include "tblite/macros.h"
#include "tblite/error.h"

/*
 * Molecular structure data class
**/

/// Molecular structure data class
typedef struct _tblite_structure* tblite_structure;

/// Create new molecular structure data (quantities in Bohr)
TBLITE_API_ENTRY tblite_structure TBLITE_API_CALL
tblite_new_structure(tblite_error /* err */,
                     const int /* natoms */,
                     const int* /* numbers [natoms] */,
                     const double* /* positions [natoms][3] */,
                     const double* /* charge in e */,
                     const int* /* uhf */,
                     const double* /* lattice [3][3] */,
                     const bool* /* periodic [3] */);

/// Delete molecular structure data
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_delete_structure(tblite_structure* /* mol */);

/// Update coordinates and lattice parameters (quantities in Bohr)
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_update_structure_geometry(tblite_error /* error */,
                                 tblite_structure /* mol */,
                                 const double* /* positions [natoms][3] */,
                                 const double* /* lattice [3][3] */);

/// Update total charge in structure object
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_update_structure_charge(tblite_error /* error */,
                               tblite_structure /* mol */,
                               const double* /* charge in e */);

/// Update number of unpaired electrons in structure object
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_update_structure_uhf(tblite_error /* error */,
                            tblite_structure /* mol */,
                            const int* /* uhf */);
